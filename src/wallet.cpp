// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "wallet.h"

#include "base58.h"
#include "clamspeech.h"
#include "coincontrol.h"
#include "kernel.h"
#include "net.h"
#include "timedata.h"
#include "txdb.h"
#include "ui_interface.h"
#include "walletdb.h"
#include "util.h"


#include <boost/algorithm/string/replace.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

// Settings
int64_t nTransactionFee = MIN_TX_FEE;
int64_t nReserveBalance = 0;
int64_t nMinimumInputValue = 0;

// limit the number of small outputs we spend at once, unless we have no option
int     nMaxOutputsToSpend = 500; // each output we spend adds 147 or 148 bytes per compressed address, and 179 or 180 bytes per uncompressed address
                                  // 500 * 180 bytes = 90k bytes, which is close to the transaction size limit

extern int64_t nMaxStakeValue;
extern int64_t nSplitSize;
extern int64_t nCombineLimit;
extern bool fCombineAny;

static unsigned int GetStakeSplitAge() { return 1 * 24 * 60 * 60; }

extern vector<CKeyID> vChangeAddresses;
extern CKeyID staketokeyID, rewardtokeyID;
extern bool fStakeTo, fRewardTo;
extern set<CBitcoinAddress> setSpendLastAddresses;
extern set<CBitcoinAddress> setStakeAddresses;

extern map<string, CClamour*> mapClamour;

int64_t gcd(int64_t n,int64_t m) { return m == 0 ? n : gcd(m, n % m); }
static inline uint64_t CoinWeightCost(const COutput &out)
{
    int64_t nTimeWeight = min((int64_t)GetTime() - (int64_t)out.tx->nTime, (int64_t)nStakeMaxAge);
    CBigNum bnCoinDayWeight = CBigNum(out.tx->vout[out.i].nValue) * nTimeWeight / (24 * 60 * 60);
    return bnCoinDayWeight.getuint64();
}

//////////////////////////////////////////////////////////////////////////////
//
// mapWallet
//

struct CompareValueOnly
{
    bool operator()(const pair<int64_t, pair<const CWalletTx*, unsigned int> >& t1,
                    const pair<int64_t, pair<const CWalletTx*, unsigned int> >& t2) const
    {
        return t1.first < t2.first;
    }
};

CPubKey CWallet::GenerateNewKey()
{
    AssertLockHeld(cs_wallet); // mapKeyMetadata
    bool fCompressed = CanSupportFeature(FEATURE_COMPRPUBKEY); // default to compressed public keys if we want 0.6.0 wallets

    RandAddSeedPerfmon();
    CKey secret;
    secret.MakeNewKey(fCompressed);

    // Compressed public keys were introduced in version 0.6.0
    if (fCompressed)
        SetMinVersion(FEATURE_COMPRPUBKEY);

    CPubKey pubkey = secret.GetPubKey();

    // Create new metadata
    int64_t nCreationTime = GetTime();
    mapKeyMetadata[pubkey.GetID()] = CKeyMetadata(nCreationTime);
    if (!nTimeFirstKey || nCreationTime < nTimeFirstKey)
        nTimeFirstKey = nCreationTime;

    if (!AddKeyPubKey(secret, pubkey))
        throw std::runtime_error("CWallet::GenerateNewKey() : AddKey failed");
    return pubkey;
}

bool CWallet::AddKeyPubKey(const CKey& secret, const CPubKey &pubkey)
{
    AssertLockHeld(cs_wallet); // mapKeyMetadata
    if (!CCryptoKeyStore::AddKeyPubKey(secret, pubkey))
        return false;
    if (!fFileBacked) 
        return true;
    if (!IsCrypted()) 
        return CWalletDB(strWalletFile).WriteKey(pubkey, secret.GetPrivKey(), mapKeyMetadata[pubkey.GetID()]);
    return true;
}

bool CWallet::RemovePubKey(const CPubKey &pubkey)
{
    AssertLockHeld(cs_wallet); // mapKeyMetadata
    if (!CCryptoKeyStore::RemovePubKey(pubkey))
        return false;
    if (!fFileBacked) 
        return true;
    return CWalletDB(strWalletFile).EraseKey(pubkey);
}

bool CWallet::AddCryptedKey(const CPubKey &vchPubKey, const vector<unsigned char> &vchCryptedSecret)
{
    if (!CCryptoKeyStore::AddCryptedKey(vchPubKey, vchCryptedSecret))
        return false;
    if (!fFileBacked)
        return true;
    {
        LOCK(cs_wallet);
        if (pwalletdbEncryption)
            return pwalletdbEncryption->WriteCryptedKey(vchPubKey, vchCryptedSecret, mapKeyMetadata[vchPubKey.GetID()]);
        else
            return CWalletDB(strWalletFile).WriteCryptedKey(vchPubKey, vchCryptedSecret, mapKeyMetadata[vchPubKey.GetID()]);
    }
    return false;
}

bool CWallet::LoadKeyMetadata(const CPubKey &pubkey, const CKeyMetadata &meta)
{
    AssertLockHeld(cs_wallet); // mapKeyMetadata
    if (meta.nCreateTime && (!nTimeFirstKey || meta.nCreateTime < nTimeFirstKey))
        nTimeFirstKey = meta.nCreateTime;

    mapKeyMetadata[pubkey.GetID()] = meta;
    return true;
}

bool CWallet::LoadCryptedKey(const CPubKey &vchPubKey, const std::vector<unsigned char> &vchCryptedSecret)
{
    return CCryptoKeyStore::AddCryptedKey(vchPubKey, vchCryptedSecret);
}

bool CWallet::AddCScript(const CScript& redeemScript)
{
    if (!CCryptoKeyStore::AddCScript(redeemScript))
        return false;
    if (!fFileBacked)
        return true;
    return CWalletDB(strWalletFile).WriteCScript(Hash160(redeemScript), redeemScript);
}

// optional setting to unlock wallet for staking only
// serves to disable the trivial sendmoney when OS account compromised
// provides no real security
bool fWalletUnlockStakingOnly = false;

bool CWallet::LoadCScript(const CScript& redeemScript)
{
    /* A sanity check was added in pull #3843 to avoid adding redeemScripts
     * that never can be redeemed. However, old wallets may still contain
     * these. Do not add them to the wallet and warn. */
    if (redeemScript.size() > MAX_SCRIPT_ELEMENT_SIZE)
    {
        std::string strAddr = CBitcoinAddress(redeemScript.GetID()).ToString();
        LogPrintf("%s: Warning: This wallet contains a redeemScript of size %"PRIszu" which exceeds maximum size %i thus can never be redeemed. Do not use address %s.\n",
            __func__, redeemScript.size(), MAX_SCRIPT_ELEMENT_SIZE, strAddr);
        return true;
    }

    return CCryptoKeyStore::AddCScript(redeemScript);
}

bool CWallet::Unlock(const SecureString& strWalletPassphrase)
{
    CCrypter crypter;
    CKeyingMaterial vMasterKey;

    {
        LOCK(cs_wallet);
        BOOST_FOREACH(const MasterKeyMap::value_type& pMasterKey, mapMasterKeys)
        {
            if(!crypter.SetKeyFromPassphrase(strWalletPassphrase, pMasterKey.second.vchSalt, pMasterKey.second.nDeriveIterations, pMasterKey.second.nDerivationMethod))
                return false;
            if (!crypter.Decrypt(pMasterKey.second.vchCryptedKey, vMasterKey))
                continue; // try another master key
            if (CCryptoKeyStore::Unlock(vMasterKey))
                return true;
        }
    }
    return false;
}

bool CWallet::ChangeWalletPassphrase(const SecureString& strOldWalletPassphrase, const SecureString& strNewWalletPassphrase)
{
    bool fWasLocked = IsLocked();

    {
        LOCK(cs_wallet);
        Lock();

        CCrypter crypter;
        CKeyingMaterial vMasterKey;
        BOOST_FOREACH(MasterKeyMap::value_type& pMasterKey, mapMasterKeys)
        {
            if(!crypter.SetKeyFromPassphrase(strOldWalletPassphrase, pMasterKey.second.vchSalt, pMasterKey.second.nDeriveIterations, pMasterKey.second.nDerivationMethod))
                return false;
            if (!crypter.Decrypt(pMasterKey.second.vchCryptedKey, vMasterKey))
                return false;
            if (CCryptoKeyStore::Unlock(vMasterKey))
            {
                int64_t nStartTime = GetTimeMillis();
                crypter.SetKeyFromPassphrase(strNewWalletPassphrase, pMasterKey.second.vchSalt, pMasterKey.second.nDeriveIterations, pMasterKey.second.nDerivationMethod);
                pMasterKey.second.nDeriveIterations = pMasterKey.second.nDeriveIterations * (100 / ((double)(GetTimeMillis() - nStartTime)));

                nStartTime = GetTimeMillis();
                crypter.SetKeyFromPassphrase(strNewWalletPassphrase, pMasterKey.second.vchSalt, pMasterKey.second.nDeriveIterations, pMasterKey.second.nDerivationMethod);
                pMasterKey.second.nDeriveIterations = (pMasterKey.second.nDeriveIterations + pMasterKey.second.nDeriveIterations * 100 / ((double)(GetTimeMillis() - nStartTime))) / 2;

                if (pMasterKey.second.nDeriveIterations < 25000)
                    pMasterKey.second.nDeriveIterations = 25000;

                LogPrintf("Wallet passphrase changed to an nDeriveIterations of %i\n", pMasterKey.second.nDeriveIterations);

                if (!crypter.SetKeyFromPassphrase(strNewWalletPassphrase, pMasterKey.second.vchSalt, pMasterKey.second.nDeriveIterations, pMasterKey.second.nDerivationMethod))
                    return false;
                if (!crypter.Encrypt(vMasterKey, pMasterKey.second.vchCryptedKey))
                    return false;
                CWalletDB(strWalletFile).WriteMasterKey(pMasterKey.first, pMasterKey.second);
                if (fWasLocked)
                    Lock();
                return true;
            }
        }
    }

    return false;
}

void CWallet::SetBestChain(const CBlockLocator& loc)
{
    CWalletDB walletdb(strWalletFile);
    walletdb.WriteBestBlock(loc);
}

bool CWallet::SetMinVersion(enum WalletFeature nVersion, CWalletDB* pwalletdbIn, bool fExplicit)
{
    LOCK(cs_wallet); // nWalletVersion
    if (nWalletVersion >= nVersion)
        return true;

    // when doing an explicit upgrade, if we pass the max version permitted, upgrade all the way
    if (fExplicit && nVersion > nWalletMaxVersion)
            nVersion = FEATURE_LATEST;

    nWalletVersion = nVersion;

    if (nVersion > nWalletMaxVersion)
        nWalletMaxVersion = nVersion;

    if (fFileBacked)
    {
        CWalletDB* pwalletdb = pwalletdbIn ? pwalletdbIn : new CWalletDB(strWalletFile);
        if (nWalletVersion > 40000)
            pwalletdb->WriteMinVersion(nWalletVersion);
        if (!pwalletdbIn)
            delete pwalletdb;
    }

    return true;
}

bool CWallet::SetMaxVersion(int nVersion)
{
    LOCK(cs_wallet); // nWalletVersion, nWalletMaxVersion
    // cannot downgrade below current version
    if (nWalletVersion > nVersion)
        return false;

    nWalletMaxVersion = nVersion;

    return true;
}

bool CWallet::EncryptWallet(const SecureString& strWalletPassphrase)
{
    if (IsCrypted())
        return false;

    CKeyingMaterial vMasterKey;
    RandAddSeedPerfmon();

    vMasterKey.resize(WALLET_CRYPTO_KEY_SIZE);
    RAND_bytes(&vMasterKey[0], WALLET_CRYPTO_KEY_SIZE);

    CMasterKey kMasterKey(nDerivationMethodIndex);

    RandAddSeedPerfmon();
    kMasterKey.vchSalt.resize(WALLET_CRYPTO_SALT_SIZE);
    RAND_bytes(&kMasterKey.vchSalt[0], WALLET_CRYPTO_SALT_SIZE);

    CCrypter crypter;
    int64_t nStartTime = GetTimeMillis();
    crypter.SetKeyFromPassphrase(strWalletPassphrase, kMasterKey.vchSalt, 25000, kMasterKey.nDerivationMethod);
    kMasterKey.nDeriveIterations = 2500000 / ((double)(GetTimeMillis() - nStartTime));

    nStartTime = GetTimeMillis();
    crypter.SetKeyFromPassphrase(strWalletPassphrase, kMasterKey.vchSalt, kMasterKey.nDeriveIterations, kMasterKey.nDerivationMethod);
    kMasterKey.nDeriveIterations = (kMasterKey.nDeriveIterations + kMasterKey.nDeriveIterations * 100 / ((double)(GetTimeMillis() - nStartTime))) / 2;

    if (kMasterKey.nDeriveIterations < 25000)
        kMasterKey.nDeriveIterations = 25000;

    LogPrintf("Encrypting Wallet with an nDeriveIterations of %i\n", kMasterKey.nDeriveIterations);

    if (!crypter.SetKeyFromPassphrase(strWalletPassphrase, kMasterKey.vchSalt, kMasterKey.nDeriveIterations, kMasterKey.nDerivationMethod))
        return false;
    if (!crypter.Encrypt(vMasterKey, kMasterKey.vchCryptedKey))
        return false;

    {
        LOCK(cs_wallet);
        mapMasterKeys[++nMasterKeyMaxID] = kMasterKey;
        if (fFileBacked)
        {
            pwalletdbEncryption = new CWalletDB(strWalletFile);
            if (!pwalletdbEncryption->TxnBegin())
                return false;
            pwalletdbEncryption->WriteMasterKey(nMasterKeyMaxID, kMasterKey);
        }

        if (!EncryptKeys(vMasterKey))
        {
            if (fFileBacked)
                pwalletdbEncryption->TxnAbort();
            exit(1); //We now probably have half of our keys encrypted in memory, and half not...die and let the user reload their unencrypted wallet.
        }

        // Encryption was introduced in version 0.4.0
        SetMinVersion(FEATURE_WALLETCRYPT, pwalletdbEncryption, true);

        if (fFileBacked)
        {
            if (!pwalletdbEncryption->TxnCommit())
                exit(1); //We now have keys encrypted in memory, but no on disk...die to avoid confusion and let the user reload their unencrypted wallet.

            delete pwalletdbEncryption;
            pwalletdbEncryption = NULL;
        }

        Lock();
        Unlock(strWalletPassphrase);
        NewKeyPool();
        Lock();

        // Need to completely rewrite the wallet file; if we don't, bdb might keep
        // bits of the unencrypted private key in slack space in the database file.
        CDB::Rewrite(strWalletFile);

    }
    NotifyStatusChanged(this);

    return true;
}

int64_t CWallet::IncOrderPosNext(CWalletDB *pwalletdb)
{
    AssertLockHeld(cs_wallet); // nOrderPosNext
    int64_t nRet = nOrderPosNext++;
    if (pwalletdb) {
        pwalletdb->WriteOrderPosNext(nOrderPosNext);
    } else {
        CWalletDB(strWalletFile).WriteOrderPosNext(nOrderPosNext);
    }
    return nRet;
}

CWallet::TxItems CWallet::OrderedTxItems(std::list<CAccountingEntry>& acentries, std::string strAccount)
{
    AssertLockHeld(cs_wallet); // mapWallet
    CWalletDB walletdb(strWalletFile);

    // First: get all CWalletTx and CAccountingEntry into a sorted-by-order multimap.
    TxItems txOrdered;

    // Note: maintaining indices in the database of (account,time) --> txid and (account, time) --> acentry
    // would make this much faster for applications that do this a lot.
    for (map<uint256, CWalletTx>::iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
    {
        CWalletTx* wtx = &((*it).second);
        txOrdered.insert(make_pair(wtx->nOrderPos, TxPair(wtx, (CAccountingEntry*)0)));
    }
    acentries.clear();
    walletdb.ListAccountCreditDebit(strAccount, acentries);
    BOOST_FOREACH(CAccountingEntry& entry, acentries)
    {
        txOrdered.insert(make_pair(entry.nOrderPos, TxPair((CWalletTx*)0, &entry)));
    }

    return txOrdered;
}

void CWallet::WalletUpdateSpent(const CTransaction &tx, bool fBlock)
{
    // Anytime a signature is successfully verified, it's proof the outpoint is spent.
    // Update the wallet spent flag if it doesn't know due to wallet.dat being
    // restored from backup or the user making copies of wallet.dat.
    {
        LOCK(cs_wallet);
        BOOST_FOREACH(const CTxIn& txin, tx.vin)
        {
            map<uint256, CWalletTx>::iterator mi = mapWallet.find(txin.prevout.hash);
            if (mi != mapWallet.end())
            {
                CWalletTx& wtx = (*mi).second;
                if (txin.prevout.n >= wtx.vout.size())
                    LogPrintf("WalletUpdateSpent: bad wtx %s\n", wtx.GetHash().ToString());
                else if (!wtx.IsSpent(txin.prevout.n) && IsMine(wtx.vout[txin.prevout.n]))
                {
                    LogPrintf("WalletUpdateSpent found spent coin %s CLAM %s\n", FormatMoney(GetCredit(wtx.vout[txin.prevout.n])), wtx.GetHash().ToString());
                    wtx.MarkSpent(txin.prevout.n);
                    wtx.WriteToDisk();
                    NotifyTransactionChanged(this, txin.prevout.hash, CT_UPDATED);
                }
            }
        }

        if (fBlock)
        {
            uint256 hash = tx.GetHash();
            map<uint256, CWalletTx>::iterator mi = mapWallet.find(hash);
            CWalletTx& wtx = (*mi).second;

            bool fMine = false;
            BOOST_FOREACH(const CTxOut& txout, tx.vout)
            {
                if (IsMine(txout))
                {
                    wtx.MarkUnspent(&txout - &tx.vout[0]);
                    fMine = true;
                }
            }

            if (fMine) {
                wtx.WriteToDisk();
                NotifyTransactionChanged(this, hash, CT_UPDATED);
            }
        }

    }
}

void CWallet::MarkDirty()
{
    {
        LOCK(cs_wallet);
        BOOST_FOREACH(PAIRTYPE(const uint256, CWalletTx)& item, mapWallet)
            item.second.MarkDirty();
    }

    fAddressRewardsReady = false;
}

bool CWallet::AddToWallet(const CWalletTx& wtxIn)
{
    uint256 hash = wtxIn.GetHash();
    {
        LOCK(cs_wallet);
        // Inserts only if not already there, returns tx inserted or tx found
        pair<map<uint256, CWalletTx>::iterator, bool> ret = mapWallet.insert(make_pair(hash, wtxIn));
        CWalletTx& wtx = (*ret.first).second;
        wtx.BindWallet(this);
        bool fInsertedNew = ret.second;
        if (fInsertedNew)
        {
            wtx.nTimeReceived = GetAdjustedTime();
            wtx.nOrderPos = IncOrderPosNext();

            wtx.nTimeSmart = wtx.nTimeReceived;
            if (wtxIn.hashBlock != 0)
            {
                if (mapBlockIndex.count(wtxIn.hashBlock))
                {
                    unsigned int latestNow = wtx.nTimeReceived;
                    unsigned int latestEntry = 0;
                    {
                        // Tolerate times up to the last timestamp in the wallet not more than 5 minutes into the future
                        int64_t latestTolerated = latestNow + 300;
                        std::list<CAccountingEntry> acentries;
                        TxItems txOrdered = OrderedTxItems(acentries);
                        for (TxItems::reverse_iterator it = txOrdered.rbegin(); it != txOrdered.rend(); ++it)
                        {
                            CWalletTx *const pwtx = (*it).second.first;
                            if (pwtx == &wtx)
                                continue;
                            CAccountingEntry *const pacentry = (*it).second.second;
                            int64_t nSmartTime;
                            if (pwtx)
                            {
                                nSmartTime = pwtx->nTimeSmart;
                                if (!nSmartTime)
                                    nSmartTime = pwtx->nTimeReceived;
                            }
                            else
                                nSmartTime = pacentry->nTime;
                            if (nSmartTime <= latestTolerated)
                            {
                                latestEntry = nSmartTime;
                                if (nSmartTime > latestNow)
                                    latestNow = nSmartTime;
                                break;
                            }
                        }
                    }

                    unsigned int& blocktime = mapBlockIndex[wtxIn.hashBlock]->nTime;
                    wtx.nTimeSmart = std::max(latestEntry, std::min(blocktime, latestNow));
                }
                else
                    LogPrintf("AddToWallet() : found %s in block %s not in index\n",
                             wtxIn.GetHash().ToString(),
                             wtxIn.hashBlock.ToString());
            }
        }

        bool fUpdated = false;
        if (!fInsertedNew)
        {
            // Merge
            if (wtxIn.hashBlock != 0 && wtxIn.hashBlock != wtx.hashBlock)
            {
                wtx.hashBlock = wtxIn.hashBlock;
                fUpdated = true;
            }
            if (wtxIn.nIndex != -1 && (wtxIn.vMerkleBranch != wtx.vMerkleBranch || wtxIn.nIndex != wtx.nIndex))
            {
                wtx.vMerkleBranch = wtxIn.vMerkleBranch;
                wtx.nIndex = wtxIn.nIndex;
                fUpdated = true;
            }
            if (wtxIn.fFromMe && wtxIn.fFromMe != wtx.fFromMe)
            {
                wtx.fFromMe = wtxIn.fFromMe;
                fUpdated = true;
            }
            fUpdated |= wtx.UpdateSpent(wtxIn.vfSpent);
        }

        //// debug print
        LogPrintf("AddToWallet %s  %s%s\n", wtxIn.GetHash().ToString(), (fInsertedNew ? "new" : ""), (fUpdated ? "update" : ""));

        // Write to disk
        if (fInsertedNew || fUpdated)
            if (!wtx.WriteToDisk())
                return false;

        if (!fHaveGUI) {
            // If default receiving address gets used, replace it with a new one
            if (vchDefaultKey.IsValid()) {
                CScript scriptDefaultKey;
                scriptDefaultKey.SetDestination(vchDefaultKey.GetID());
                BOOST_FOREACH(const CTxOut& txout, wtx.vout)
                {
                    if (txout.scriptPubKey == scriptDefaultKey)
                    {
                        CPubKey newDefaultKey;
                        if (GetKeyFromPool(newDefaultKey))
                        {
                            SetDefaultKey(newDefaultKey);
                            SetAddressBookName(vchDefaultKey.GetID(), "");
                        }
                    }
                }
            }
        }
        // since AddToWallet is called directly for self-originating transactions, check for consumption of own coins
        WalletUpdateSpent(wtx, (wtxIn.hashBlock != 0));

        // Notify UI of new or updated transaction
        NotifyTransactionChanged(this, hash, fInsertedNew ? CT_NEW : CT_UPDATED);

        // notify an external script when a wallet transaction comes in or is updated
        std::string strCmd = GetArg("-walletnotify", "");

        if ( !strCmd.empty())
        {
            boost::replace_all(strCmd, "%s", wtxIn.GetHash().GetHex());
            boost::thread t(runCommand, strCmd); // thread runs free
        }

    }
    return true;
}

// Add a transaction to the wallet, or update it.
// pblock is optional, but should be provided if the transaction is known to be in a block.
// If fUpdate is true, existing transactions will be updated.
bool CWallet::AddToWalletIfInvolvingMe(const CTransaction& tx, const CBlock* pblock, bool fUpdate)
{
    uint256 hash = tx.GetHash();
    {
        LOCK(cs_wallet);
        bool fExisted = mapWallet.count(hash);
        if (fExisted && !fUpdate) return false;
        if (fExisted || IsMine(tx) || IsFromMe(tx))
        {
            CWalletTx wtx(this,tx);
            // Get merkle branch if transaction was found in a block
            if (pblock)
                wtx.SetMerkleBranch(pblock);
            return AddToWallet(wtx);
        }
        else
            WalletUpdateSpent(tx);
    }
    return false;
}

void CWallet::SyncTransaction(const CTransaction& tx, const CBlock* pblock, bool fConnect) {
    if (!fConnect)
    {
        // wallets need to refund inputs when disconnecting coinstake
        if (tx.IsCoinStake())
        {
            if (IsFromMe(tx))
                DisableTransaction(tx);
        }
        return;
    }

    AddToWalletIfInvolvingMe(tx, pblock, true);
}

void CWallet::StakeTransaction(const CScript& script, int64_t nStakeReward, bool fConnect) {
    if (!fAddressRewardsReady)
        return;

    // if (!fConnect)
    //     LogPrintf("stake tx disconnecting %s\n", FormatMoney(nStakeReward));
    // else
    //     LogPrintf("stake tx connecting    %s\n", FormatMoney(nStakeReward));

    CTxDestination address;
    if (!ExtractDestination(script, address)) {
        LogPrintf("can't extract destination\n");
        return;
    }

    std::string addr(CBitcoinAddress(address).ToString());
    if (!::IsMine(*this, address)) {
        LogPrintf("stake %s for %s; not mine\n", FormatMoney(nStakeReward), addr);
        return;
    }

    mapAddressRewards["*"] += nStakeReward;
    mapAddressRewards[addr] += nStakeReward;
    LogPrintf("stake %s for %s; global = %s, address = %s\n",
              FormatMoney(nStakeReward), addr,
              FormatMoney(mapAddressRewards["*"]),
              FormatMoney(mapAddressRewards[addr]));

    std::string strCmd = GetArg("-stakenotify", "");

    if (!strCmd.empty())
    {
        boost::replace_all(strCmd, "%r", FormatMoney(nStakeReward));
        boost::replace_all(strCmd, "%a", addr);
        boost::replace_all(strCmd, "%t", FormatMoney(mapAddressRewards["*"]));
        boost::replace_all(strCmd, "%s", FormatMoney(mapAddressRewards[addr]));

        boost::thread t(runCommand, strCmd); // thread runs free
    }
}

void CWallet::EraseFromWallet(const uint256 &hash)
{
    if (!fFileBacked)
        return;
    {
        LOCK(cs_wallet);
        if (mapWallet.erase(hash))
            CWalletDB(strWalletFile).EraseTx(hash);
    }
    return;
}


bool CWallet::IsMine(const CTxIn &txin) const
{
    {
        LOCK(cs_wallet);
        map<uint256, CWalletTx>::const_iterator mi = mapWallet.find(txin.prevout.hash);
        if (mi != mapWallet.end())
        {
            const CWalletTx& prev = (*mi).second;
            if (txin.prevout.n < prev.vout.size())
                if (IsMine(prev.vout[txin.prevout.n]))
                    return true;
        }
    }
    return false;
}

int64_t CWallet::GetDebit(const CTxIn &txin) const
{
    {
        LOCK(cs_wallet);
        map<uint256, CWalletTx>::const_iterator mi = mapWallet.find(txin.prevout.hash);
        if (mi != mapWallet.end())
        {
            const CWalletTx& prev = (*mi).second;
            if (txin.prevout.n < prev.vout.size())
                if (IsMine(prev.vout[txin.prevout.n]))
                    return prev.vout[txin.prevout.n].nValue;
        }
    }
    return 0;
}

bool CWallet::IsChange(const CTxOut& txout) const
{
    CTxDestination address;

    // TODO: fix handling of 'change' outputs. The assumption is that any
    // payment to a TX_PUBKEYHASH that is mine but isn't in the address book
    // is change. That assumption is likely to break when we implement multisignature
    // wallets that return change back into a multi-signature-protected address;
    // a better way of identifying which outputs are 'the send' and which are
    // 'the change' will need to be implemented (maybe extend CWalletTx to remember
    // which output, if any, was change).
    if (ExtractDestination(txout.scriptPubKey, address) && ::IsMine(*this, address))
    {
        LOCK(cs_wallet);
        if (!mapAddressBook.count(address))
            return true;
    }
    return false;
}

int64_t CWalletTx::GetTxTime() const
{
    int64_t n = nTimeSmart;
    return n ? n : nTimeReceived;
}

int CWalletTx::GetRequestCount() const
{
    // Returns -1 if it wasn't being tracked
    int nRequests = -1;
    {
        LOCK(pwallet->cs_wallet);
        if (IsCoinBase() || IsCoinStake())
        {
            // Generated block
            if (hashBlock != 0)
            {
                map<uint256, int>::const_iterator mi = pwallet->mapRequestCount.find(hashBlock);
                if (mi != pwallet->mapRequestCount.end())
                    nRequests = (*mi).second;
            }
        }
        else
        {
            // Did anyone request this transaction?
            map<uint256, int>::const_iterator mi = pwallet->mapRequestCount.find(GetHash());
            if (mi != pwallet->mapRequestCount.end())
            {
                nRequests = (*mi).second;

                // How about the block it's in?
                if (nRequests == 0 && hashBlock != 0)
                {
                    map<uint256, int>::const_iterator mi = pwallet->mapRequestCount.find(hashBlock);
                    if (mi != pwallet->mapRequestCount.end())
                        nRequests = (*mi).second;
                    else
                        nRequests = 1; // If it's in someone else's block it must have got out
                }
            }
        }
    }
    return nRequests;
}

void CWalletTx::GetAmounts(list<pair<CTxDestination, int64_t> >& listReceived,
                           list<pair<CTxDestination, int64_t> >& listSent, int64_t& nFee, string& strSentAccount) const
{
    nFee = 0;
    listReceived.clear();
    listSent.clear();
    strSentAccount = strFromAccount;

    // Compute fee:
    int64_t nDebit = GetDebit();
    if (nDebit > 0) // debit>0 means we signed/sent this transaction
    {
        int64_t nValueOut = GetValueOut();
        nFee = nDebit - nValueOut;
    }

    // Sent/received.
    BOOST_FOREACH(const CTxOut& txout, vout)
    {
        // Skip special stake out
        if (txout.scriptPubKey.empty())
            continue;

        bool fIsMine;
        // Only need to handle txouts if AT LEAST one of these is true:
        //   1) they debit from us (sent)
        //   2) the output is to us (received)
        if (nDebit > 0)
        {
            // Don't report 'change' txouts
            if (pwallet->IsChange(txout))
                continue;
            fIsMine = pwallet->IsMine(txout);
        }
        else if (!(fIsMine = pwallet->IsMine(txout)))
            continue;

        // In either case, we need to get the destination address
        CTxDestination address;
        if (!ExtractDestination(txout.scriptPubKey, address))
        {
            LogPrintf("CWalletTx::GetAmounts: Unknown transaction type found, txid %s\n",
                     this->GetHash().ToString());
            address = CNoDestination();
        }

        // If we are debited by the transaction, add the output as a "sent" entry
        if (nDebit > 0)
            listSent.push_back(make_pair(address, txout.nValue));

        // If we are receiving the output, add it as a "received" entry
        if (fIsMine)
            listReceived.push_back(make_pair(address, txout.nValue));
    }

}

void CWalletTx::GetAccountAmounts(const string& strAccount, int64_t& nReceived,
                                  int64_t& nSent, int64_t& nFee) const
{
    nReceived = nSent = nFee = 0;

    int64_t allFee;
    string strSentAccount;
    list<pair<CTxDestination, int64_t> > listReceived;
    list<pair<CTxDestination, int64_t> > listSent;
    GetAmounts(listReceived, listSent, allFee, strSentAccount);

    if (strAccount == strSentAccount)
    {
        BOOST_FOREACH(const PAIRTYPE(CTxDestination,int64_t)& s, listSent)
            nSent += s.second;
        nFee = allFee;
    }
    {
        LOCK(pwallet->cs_wallet);
        BOOST_FOREACH(const PAIRTYPE(CTxDestination,int64_t)& r, listReceived)
        {
            if (pwallet->mapAddressBook.count(r.first))
            {
                map<CTxDestination, string>::const_iterator mi = pwallet->mapAddressBook.find(r.first);
                if (mi != pwallet->mapAddressBook.end() && (*mi).second == strAccount)
                    nReceived += r.second;
            }
            else if (strAccount.empty())
            {
                nReceived += r.second;
            }
        }
    }
}

void CWalletTx::AddSupportingTransactions()
{
    vtxPrev.clear();

    const int COPY_DEPTH = 3;
    if (SetMerkleBranch() < COPY_DEPTH)
    {
        vector<uint256> vWorkQueue;
        BOOST_FOREACH(const CTxIn& txin, vin)
            vWorkQueue.push_back(txin.prevout.hash);

        // This critsect is OK because txdb is already open
        {
            LOCK(pwallet->cs_wallet);
            map<uint256, const CMerkleTx*> mapWalletPrev;
            set<uint256> setAlreadyDone;
            for (unsigned int i = 0; i < vWorkQueue.size(); i++)
            {
                uint256 hash = vWorkQueue[i];
                if (setAlreadyDone.count(hash))
                    continue;
                setAlreadyDone.insert(hash);

                CMerkleTx tx;
                map<uint256, CWalletTx>::const_iterator mi = pwallet->mapWallet.find(hash);
                if (mi != pwallet->mapWallet.end())
                {
                    tx = (*mi).second;
                    BOOST_FOREACH(const CMerkleTx& txWalletPrev, (*mi).second.vtxPrev)
                        mapWalletPrev[txWalletPrev.GetHash()] = &txWalletPrev;
                }
                else if (mapWalletPrev.count(hash))
                {
                    tx = *mapWalletPrev[hash];
                }
                else
                {
                    LogPrintf("ERROR: AddSupportingTransactions() : unsupported transaction\n");
                    continue;
                }

                int nDepth = tx.SetMerkleBranch();
                vtxPrev.push_back(tx);

                if (nDepth < COPY_DEPTH)
                {
                    BOOST_FOREACH(const CTxIn& txin, tx.vin)
                        vWorkQueue.push_back(txin.prevout.hash);
                }
            }
        }
    }

    reverse(vtxPrev.begin(), vtxPrev.end());
}

bool CWalletTx::WriteToDisk()
{
    return CWalletDB(pwallet->strWalletFile).WriteTx(GetHash(), *this);
}

// Scan the block chain (starting in pindexStart) for transactions
// from or to us. If fUpdate is true, found transactions that already
// exist in the wallet will be updated.
int CWallet::ScanForWalletTransactions(CBlockIndex* pindexStart, bool fUpdate, bool fFullScan)
{
    int ret = 0;

    CBlockIndex* pindex = pindexStart;
    {
        

        LOCK2(cs_main, cs_wallet);
        while (pindex)
        {
            // no need to read and scan block, if block was created before
            // our wallet birthday (as adjusted for block time variability)
                if (nTimeFirstKey && (pindex->nTime < (nTimeFirstKey - 7200)) && !fFullScan) {
                    pindex = pindex->pnext;
                    continue;
                }
            CBlock block;
            block.ReadFromDisk(pindex);
            BOOST_FOREACH(CTransaction& tx, block.vtx)
            {
                if (AddToWalletIfInvolvingMe(tx, &block, fUpdate))
                    ret++;
            }
            pindex = pindex->pnext;
        }
    }
    return ret;
}

void CWallet::ReacceptWalletTransactions()
{
    bool fRepeat = true;
    while (fRepeat)
    {
        LOCK2(cs_main, cs_wallet);
        fRepeat = false;
        bool fMissing = false;
        BOOST_FOREACH(PAIRTYPE(const uint256, CWalletTx)& item, mapWallet)
        {
            CWalletTx& wtx = item.second;
            if ((wtx.IsCoinBase() && wtx.IsSpent(0)) || (wtx.IsCoinStake() && wtx.IsSpent(1)))
                continue;

            CCoins coins;
            bool fUpdated = false;
            bool fFound = pcoinsTip->GetCoins(wtx.GetHash(), coins);
            if (fFound || wtx.GetDepthInMainChain() > 0)
            {
                // Update fSpent if a tx got spent somewhere else by a copy of wallet.dat
                for (unsigned int i = 0; i < wtx.vout.size(); i++)
                {
                    if (wtx.IsSpent(i))
                        continue;
                    if ((i >= coins.vout.size() || coins.vout[i].IsNull()) && IsMine(wtx.vout[i]))
                    {
                        wtx.MarkSpent(i);
                        fUpdated = true;
                        fMissing = true;
                    }
                }
                if (fUpdated)
                {
                    LogPrintf("ReacceptWalletTransactions found spent coin %s CLAM %s\n", FormatMoney(wtx.GetCredit()), wtx.GetHash().ToString());
                    wtx.MarkDirty();
                    wtx.WriteToDisk();
                }
            }
            else
            {
                // Re-accept any txes of ours that aren't already in a block
                if (!(wtx.IsCoinBase() || wtx.IsCoinStake()))
                    wtx.AcceptWalletTransaction(false);
            }
        }
        if (fMissing)
        {
            // TODO: optimize this to scan just part of the block chain?
            if (ScanForWalletTransactions(pindexGenesisBlock))
                fRepeat = true;  // Found missing transactions: re-do re-accept.
        }
    }
}

void CWalletTx::RelayWalletTransaction()
{
    BOOST_FOREACH(const CMerkleTx& tx, vtxPrev)
    {
        // Important: versions of bitcoin before 0.8.6 had a bug that inserted
        // empty transactions into the vtxPrev, which will cause the node to be
        // banned when retransmitted, hence the check for !tx.vin.empty()
        if (!tx.IsCoinBase() && !tx.IsCoinStake() && !tx.vin.empty())
            if (tx.GetDepthInMainChain() == 0)
                RelayTransaction((CTransaction)tx, tx.GetHash());
    }
    if (!IsCoinBase() && !IsCoinStake())
    {
        if (GetDepthInMainChain() == 0) {
            uint256 hash = GetHash();
            LogPrintf("Relaying wtx %s\n", hash.ToString().c_str());
            RelayTransaction((CTransaction)*this, hash);
        }
    }
}

void CWallet::ResendWalletTransactions(bool fForce)
{
    // Do this infrequently and randomly to avoid giving away
    // that these are our transactions.
    static int64 nNextTime;
    if (GetTime() < nNextTime)
        return;
    bool fFirst = (nNextTime == 0);
    nNextTime = GetTime() + GetRand(30 * 60);
    if (fFirst)
        return;

    // Only do it if there's been a new block since last time
    static int64 nLastTime;
    if (nTimeBestReceived < nLastTime)
        return;
    nLastTime = GetTime();

    // Rebroadcast any of our txes that aren't in a block yet
    LogPrintf("ResendWalletTransactions()\n");
    {
        LOCK(cs_wallet);
        // Sort them in chronological order
        multimap<unsigned int, CWalletTx*> mapSorted;
        BOOST_FOREACH(PAIRTYPE(const uint256, CWalletTx)& item, mapWallet)
        {
            CWalletTx& wtx = item.second;
            // Don't rebroadcast until it's had plenty of time that
            // it should have gotten in already by now.
            if (nTimeBestReceived - (int64)wtx.nTimeReceived > 5 * 60)
                mapSorted.insert(make_pair(wtx.nTimeReceived, &wtx));
        }
        BOOST_FOREACH(PAIRTYPE(const unsigned int, CWalletTx*)& item, mapSorted)
        {
            CWalletTx& wtx = *item.second;
            wtx.RelayWalletTransaction();
        }
    }
}

void CWallet::TidyWalletTransactions()
{
    vector<uint256> vHashes;
    for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
    {
        const CWalletTx* pcoin = &(*it).second;
        if (!IsMine(*pcoin) && !IsFromMe(*pcoin))
            vHashes.push_back(pcoin->GetHash());
    }

    BOOST_FOREACH(uint256& hash, vHashes)
    {
        EraseFromWallet(hash);
        NotifyTransactionChanged(this, hash, CT_DELETED);
    }
}





//////////////////////////////////////////////////////////////////////////////
//
// Actions
//


int64_t CWallet::GetBalance() const
{
    int64_t nTotal = 0;
    {
        LOCK2(cs_main, cs_wallet);
        for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
        {
            const CWalletTx* pcoin = &(*it).second;
            if (pcoin->IsTrusted())
                nTotal += pcoin->GetAvailableCredit();
        }
    }

    return nTotal;
}

int64_t CWallet::GetUnconfirmedBalance() const
{
    int64_t nTotal = 0;
    {
        LOCK2(cs_main, cs_wallet);
        for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
        {
            const CWalletTx* pcoin = &(*it).second;
            if (!IsFinalTx(*pcoin) || (!pcoin->IsTrusted() && pcoin->GetDepthInMainChain() == 0))
                nTotal += pcoin->GetAvailableCredit();
        }
    }
    return nTotal;
}

int64_t CWallet::GetImmatureBalance() const
{
    int64_t nTotal = 0;
    {
        LOCK2(cs_main, cs_wallet);
        for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
        {
            const CWalletTx& pcoin = (*it).second;
            if (pcoin.IsCoinBase() && pcoin.GetBlocksToMaturity() > 0 && pcoin.IsInMainChain())
                nTotal += GetCredit(pcoin);
        }
    }
    return nTotal;
}

// populate vCoins with vector of spendable COutputs
void CWallet::AvailableCoins(vector<COutput>& vCoins, bool fOnlyConfirmed, const CCoinControl *coinControl, bool fOnlyMature) const
{
    vCoins.clear();

    {
        LOCK2(cs_main, cs_wallet);
        for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
        {
            const CWalletTx* pcoin = &(*it).second;

            if (!IsFinalTx(*pcoin))
                continue;

            if (fOnlyConfirmed && !pcoin->IsTrusted())
                continue;

            if (pcoin->IsCoinBase() && fOnlyMature && pcoin->GetBlocksToMaturity() > 0)
                continue;

            if(pcoin->IsCoinStake() && fOnlyMature && pcoin->GetBlocksToMaturity() > 0)
                continue;

            int nDepth = pcoin->GetDepthInMainChain();
            if (nDepth < 0)
                continue;

            for (unsigned int i = 0; i < pcoin->vout.size(); i++)
                if (!(pcoin->IsSpent(i)) && IsMine(pcoin->vout[i]) && pcoin->vout[i].nValue >= nMinimumInputValue &&
                (!coinControl || !coinControl->HasSelected() || coinControl->IsSelected((*it).first, i)))
                    vCoins.push_back(COutput(pcoin, i, nDepth));

        }
    }
}

void CWallet::AvailableCoinsForStaking(vector<COutput>& vCoins, unsigned int nSpendTime) const
{
    vCoins.clear();

    {
        LOCK2(cs_main, cs_wallet);
        for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
        {
            const CWalletTx* pcoin = &(*it).second;

            // Filtering by tx timestamp instead of block timestamp may give false positives but never false negatives
            if (pcoin->nTime + nStakeMinAge > nSpendTime)
                continue;
    
            if (pcoin->GetBlocksToMaturity() > 0)
                continue;

            int nDepth = pcoin->GetDepthInMainChain();
            if (nDepth < 1)
                continue;

            CTxDestination address;
            for (unsigned int i = 0; i < pcoin->vout.size(); i++)
                if (!(pcoin->IsSpent(i)) && IsMine(pcoin->vout[i]) && pcoin->vout[i].nValue >= nMinimumInputValue &&
                    (!nMaxStakeValue || pcoin->vout[i].nValue <= nMaxStakeValue) &&
                    (!setStakeAddresses.size() || (ExtractDestination(pcoin->vout[i].scriptPubKey, address) && setStakeAddresses.count(address)))) {
                    pcoin->hash = it->first;
                    vCoins.push_back(COutput(pcoin, i, nDepth));
                }
        }
    }
}

static void ApproximateBestSubset(vector<pair<int64_t, pair<const CWalletTx*,unsigned int> > >vValue, int64_t nTotalLower, int64_t nTargetValue,
                                  vector<char>& vfBest, int64_t& nBest, int& nSize, int iterations = 1000)
{
    vector<char> vfIncluded;

    vfBest.assign(vValue.size(), true);
    nBest = nTotalLower;

    seed_insecure_rand();

    for (int nRep = 0; nRep < iterations && nBest != nTargetValue; nRep++)
    {
        vfIncluded.assign(vValue.size(), false);
        int64_t nTotal = 0;
        bool fReachedTarget = false;
        for (int nPass = 0; nPass < 2 && !fReachedTarget; nPass++)
        {
            for (unsigned int i = 0; i < vValue.size(); i++)
            {
                //The solver here uses a randomized algorithm,
                //the randomness serves no real security purpose but is just
                //needed to prevent degenerate behavior and it is important
                //that the rng fast. We do not use a constant random sequence,
                //because there may be some privacy improvement by making
                //the selection random.
                if (nPass == 0 ? insecure_rand()&1 : !vfIncluded[i])
                {
                    nTotal += vValue[i].first;
                    vfIncluded[i] = true;
                    if (nTotal >= nTargetValue)
                    {
                        fReachedTarget = true;
                        if (nTotal < nBest)
                        {
                            nBest = nTotal;
                            vfBest = vfIncluded;
                        }
                        nTotal -= vValue[i].first;
                        vfIncluded[i] = false;
                    }
                }
            }
        }
    }

    // count the included outputs
    nSize = 0;
    for (unsigned int i = 0; i < vValue.size(); i++)
        if (vfBest[i])
            nSize++;
}

// ppcoin: total coins staked (non-spendable until maturity)
int64_t CWallet::GetStake() const
{
    int64_t nTotal = 0;
    LOCK2(cs_main, cs_wallet);
    for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
    {
        const CWalletTx* pcoin = &(*it).second;
        if (pcoin->IsCoinStake() && pcoin->GetBlocksToMaturity() > 0 && pcoin->GetDepthInMainChain() > 0)
            nTotal += CWallet::GetCredit(*pcoin);
    }
    return nTotal;
}

int64_t CWallet::GetNewMint() const
{
    int64_t nTotal = 0;
    LOCK2(cs_main, cs_wallet);
    for (map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
    {
        const CWalletTx* pcoin = &(*it).second;
        if (pcoin->IsCoinBase() && pcoin->GetBlocksToMaturity() > 0 && pcoin->GetDepthInMainChain() > 0)
            nTotal += CWallet::GetCredit(*pcoin);
    }
    return nTotal;
}

struct LargerOrEqualThanThreshold
{
    int64_t threshold;
    LargerOrEqualThanThreshold(int64_t threshold) : threshold(threshold) {}
    bool operator()(pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > const &v) const { return v.first.first >= threshold; }
};

bool CWallet::SelectCoinsMinConfByCoinAge(int64_t nTargetValue, unsigned int nSpendTime, int nConfMine, int nConfTheirs, std::vector<COutput> vCoins, set<pair<const CWalletTx*,unsigned int> >& setCoinsRet, int64_t& nValueRet, bool fLogFailureReason) const
{
    setCoinsRet.clear();
    nValueRet = 0;

    vector<pair<COutput, uint64_t> > mCoins;
    BOOST_FOREACH(const COutput& out, vCoins)
    {
        mCoins.push_back(std::make_pair(out, CoinWeightCost(out)));
    }

    // List of values less than target
    pair<pair<int64_t,int64_t>, pair<const CWalletTx*,unsigned int> > coinLowestLarger;
    coinLowestLarger.first.second = std::numeric_limits<int64_t>::max();
    coinLowestLarger.second.first = NULL;
    vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > > vValue;
    int64_t nTotalLower = 0;
    boost::sort(mCoins, boost::bind(&std::pair<COutput, uint64_t>::second, _1) < boost::bind(&std::pair<COutput, uint64_t>::second, _2));

    BOOST_FOREACH(const PAIRTYPE(COutput, uint64_t)& output, mCoins)
    {
        const CWalletTx *pcoin = output.first.tx;

        if (output.first.nDepth < (pcoin->IsFromMe() ? nConfMine : nConfTheirs))
            continue;

        int i = output.first.i;

        // Follow the timestamp rules
        if (pcoin->nTime > nSpendTime)
            continue;

        int64_t n = pcoin->vout[i].nValue;

        pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > coin = make_pair(make_pair(n,output.second),make_pair(pcoin, i));

        if (n < nTargetValue + CENT)
        {
            vValue.push_back(coin);
            nTotalLower += n;
        }
        else if (output.second < (uint64_t)coinLowestLarger.first.second)
        {
            coinLowestLarger = coin;
        }
    }

    if (nTotalLower < nTargetValue)
    {
        if (coinLowestLarger.second.first == NULL) {
            if (fLogFailureReason)
                LogPrintf("CWallet::SelectCoinsMinConfByCoinAge failed: sum of lower outputs (%s) < target (%s), and no other outputs\n",
                          FormatMoney(nTotalLower), FormatMoney(nTargetValue));
            return false;
        }
        setCoinsRet.insert(coinLowestLarger.second);
        nValueRet += coinLowestLarger.first.first;
        return true;
    }

    // Calculate dynamic programming matrix
    int64_t nTotalValue = vValue[0].first.first;
    int64_t nGCD = vValue[0].first.first;
    for (unsigned int i = 1; i < vValue.size(); ++i)
    {
        nGCD = gcd(vValue[i].first.first, nGCD);
        nTotalValue += vValue[i].first.first;
    }
    nGCD = gcd(nTargetValue, nGCD);
    int64_t denom = nGCD;
    const int64_t k = 25;
    const int64_t approx = int64_t(vValue.size() * (nTotalValue - nTargetValue)) / k;
    if (approx > nGCD)
    {
        denom = approx; // apply approximation
    }
    if (fDebug) cerr << "nGCD " << nGCD << " denom " << denom << " k " << k << endl;

    if (nTotalValue == nTargetValue)
    {
        for (unsigned int i = 0; i < vValue.size(); ++i)
        {
            setCoinsRet.insert(vValue[i].second);
        }
        nValueRet = nTotalValue;
        return true;
    }

    size_t nBeginBundles = vValue.size();
    size_t nTotalCoinValues = vValue.size();
    size_t nBeginCoinValues = 0;
    int64_t costsum = 0;
    vector<vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > >::iterator> vZeroValueBundles;
    if (denom != nGCD)
    {
        // All coin outputs that with zero value will always be added by the dynamic programming routine
        // So we collect them into bundles of value denom
        vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > >::iterator itZeroValue = std::stable_partition(vValue.begin(), vValue.end(), LargerOrEqualThanThreshold(denom));
        vZeroValueBundles.push_back(itZeroValue);
        pair<int64_t, int64_t> pBundle = make_pair(0, 0);
        nBeginBundles = itZeroValue - vValue.begin();
        nTotalCoinValues = nBeginBundles;
        while (itZeroValue != vValue.end())
        {
            pBundle.first += itZeroValue->first.first;
            pBundle.second += itZeroValue->first.second;
            itZeroValue++;
            if (pBundle.first >= denom)
            {
                vZeroValueBundles.push_back(itZeroValue);
                vValue[nTotalCoinValues].first = pBundle;
                pBundle = make_pair(0, 0);
                nTotalCoinValues++;
            }
        }
        // We need to recalculate the total coin value due to truncation of integer division
        nTotalValue = 0;
        for (unsigned int i = 0; i < nTotalCoinValues; ++i)
        {
            nTotalValue += vValue[i].first.first / denom;
        }
        // Check if dynamic programming is still applicable with the approximation
        if (nTargetValue/denom >= nTotalValue)
        {
            // We lose too much coin value through the approximation, i.e. the residual of the previous recalculation is too large
            // Since the partitioning of the previously sorted list is stable, we can just pick the first coin outputs in the list until we have a valid target value
            for (; nBeginCoinValues < nTotalCoinValues && (nTargetValue - nValueRet)/denom >= nTotalValue; ++nBeginCoinValues)
            {
                if (nBeginCoinValues >= nBeginBundles)
                {
                    if (fDebug) cerr << "prepick bundle item " << FormatMoney(vValue[nBeginCoinValues].first.first) << " normalized " << vValue[nBeginCoinValues].first.first / denom << " cost " << vValue[nBeginCoinValues].first.second << endl;
                    const size_t nBundle = nBeginCoinValues - nBeginBundles;
                    for (vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > >::iterator it = vZeroValueBundles[nBundle]; it != vZeroValueBundles[nBundle + 1]; ++it)
                    {
                        setCoinsRet.insert(it->second);
                    }
                }
                else
                {
                    if (fDebug) cerr << "prepicking " << FormatMoney(vValue[nBeginCoinValues].first.first) << " normalized " << vValue[nBeginCoinValues].first.first / denom << " cost " << vValue[nBeginCoinValues].first.second << endl;
                    setCoinsRet.insert(vValue[nBeginCoinValues].second);
                }
                nTotalValue -= vValue[nBeginCoinValues].first.first / denom;
                nValueRet += vValue[nBeginCoinValues].first.first;
                costsum += vValue[nBeginCoinValues].first.second;
            }
            if (nValueRet >= nTargetValue)
            {
                    if (fDebug) cerr << "Done without dynprog: " << "requested " << FormatMoney(nTargetValue) << "\tnormalized " << nTargetValue/denom + (nTargetValue % denom != 0 ? 1 : 0) << "\tgot " << FormatMoney(nValueRet) << "\tcost " << costsum << endl;
                    return true;
            }
        }
    }
    else
    {
        nTotalValue /= denom;
    }

    uint64_t nAppend = 1;
    if ((nTargetValue - nValueRet) % denom != 0)
    {
        // We need to decrease the capacity because of integer truncation
        nAppend--;
    }

    // The capacity (number of columns) corresponds to the amount of coin value we are allowed to discard
    boost::numeric::ublas::matrix<uint64_t> M((nTotalCoinValues - nBeginCoinValues) + 1, (nTotalValue - (nTargetValue - nValueRet)/denom) + nAppend, std::numeric_limits<int64_t>::max());
    boost::numeric::ublas::matrix<unsigned int> B((nTotalCoinValues - nBeginCoinValues) + 1, (nTotalValue - (nTargetValue - nValueRet)/denom) + nAppend);
    for (unsigned int j = 0; j < M.size2(); ++j)
    {
        M(0,j) = 0;
    }
    for (unsigned int i = 1; i < M.size1(); ++i)
    {
        uint64_t nWeight = vValue[nBeginCoinValues + i - 1].first.first / denom;
        uint64_t nValue = vValue[nBeginCoinValues + i - 1].first.second;
        //cerr << "Weight " << nWeight << " Value " << nValue << endl;
        for (unsigned int j = 0; j < M.size2(); ++j)
        {
            B(i, j) = j;
            if (nWeight <= j)
            {
                uint64_t nStep = M(i - 1, j - nWeight) + nValue;
                if (M(i - 1, j) >= nStep)
                {
                    M(i, j) = M(i - 1, j);
                }
                else
                {
                    M(i, j) = nStep;
                    B(i, j) = j - nWeight;
                }
            }
            else
            {
                M(i, j) = M(i - 1, j);
            }
        }
    }
    // Trace back optimal solution
    int64_t nPrev = M.size2() - 1;
    for (unsigned int i = M.size1() - 1; i > 0; --i)
    {
        //cerr << i - 1 << " " << vValue[i - 1].second.second << " " << vValue[i - 1].first.first << " " << vValue[i - 1].first.second << " " << nTargetValue << " " << nPrev << " " << (nPrev == B(i, nPrev) ? "XXXXXXXXXXXXXXX" : "") << endl;
        if (nPrev == B(i, nPrev))
        {
            const size_t nValue = nBeginCoinValues + i - 1;
            // Check if this is a bundle
            if (nValue >= nBeginBundles)
            {
                if (fDebug) cerr << "pick bundle item " << FormatMoney(vValue[nValue].first.first) << " normalized " << vValue[nValue].first.first / denom << " cost " << vValue[nValue].first.second << endl;
                const size_t nBundle = nValue - nBeginBundles;
                for (vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > >::iterator it = vZeroValueBundles[nBundle]; it != vZeroValueBundles[nBundle + 1]; ++it)
                {
                    setCoinsRet.insert(it->second);
                }
            }
            else
            {
                if (fDebug) cerr << "pick " << nValue << " value " << FormatMoney(vValue[nValue].first.first) << " normalized " << vValue[nValue].first.first / denom << " cost " << vValue[nValue].first.second << endl;
                setCoinsRet.insert(vValue[nValue].second);
            }
            nValueRet += vValue[nValue].first.first;
            costsum += vValue[nValue].first.second;
        }
        nPrev = B(i, nPrev);
    }
    if (nValueRet < nTargetValue && !vZeroValueBundles.empty())
    {
        // If we get here it means that there are either not sufficient funds to pay the transaction or that there are small coin outputs left that couldn't be bundled
        // We try to fulfill the request by adding these small coin outputs
        for (vector<pair<pair<int64_t,int64_t>,pair<const CWalletTx*,unsigned int> > >::iterator it = vZeroValueBundles.back(); it != vValue.end() && nValueRet < nTargetValue; ++it)
        {
             setCoinsRet.insert(it->second);
             nValueRet += it->first.first;
        }
    }
    if (fDebug) cerr << "requested " << FormatMoney(nTargetValue) << "\tnormalized " << nTargetValue/denom + (nTargetValue % denom != 0 ? 1 : 0) << "\tgot " << FormatMoney(nValueRet) << "\tcost " << costsum << endl;
    if (fDebug) cerr << "M " << M.size1() << "x" << M.size2() << "; vValue.size() = " << vValue.size() << endl;
    return true;
}

bool CWallet::SelectCoinsMinConf(int64_t nTargetValue, unsigned int nSpendTime, int nConfMine, int nConfTheirs, vector<COutput> vCoins, set<pair<const CWalletTx*,unsigned int> >& setCoinsRet, int64_t& nValueRet, bool fLogFailureReason) const
{
    setCoinsRet.clear();
    nValueRet = 0;

    // List of values less than target
    pair<int64_t, pair<const CWalletTx*,unsigned int> > coinLowestLarger;
    coinLowestLarger.first = std::numeric_limits<int64_t>::max();
    coinLowestLarger.second.first = NULL;
    vector<pair<int64_t, pair<const CWalletTx*,unsigned int> > > vValue;
    int64_t nTotalLower = 0;

    random_shuffle(vCoins.begin(), vCoins.end(), GetRandInt);

    BOOST_FOREACH(COutput output, vCoins)
    {
        const CWalletTx *pcoin = output.tx;

        if (output.nDepth < (pcoin->IsFromMe() ? nConfMine : nConfTheirs))
            continue;

        int i = output.i;

        // Follow the timestamp rules
        if (pcoin->nTime > nSpendTime)
            continue;

        int64_t n = pcoin->vout[i].nValue;

        pair<int64_t,pair<const CWalletTx*,unsigned int> > coin = make_pair(n,make_pair(pcoin, i));

        if (n == nTargetValue)
        {
            setCoinsRet.insert(coin.second);
            nValueRet += coin.first;
            return true;
        }
        else if (n < nTargetValue + CENT)
        {
            vValue.push_back(coin);
            nTotalLower += n;
        }
        else if (n < coinLowestLarger.first)
        {
            coinLowestLarger = coin;
        }
    }

    if (nTotalLower == nTargetValue)
    {
        for (unsigned int i = 0; i < vValue.size(); ++i)
        {
            setCoinsRet.insert(vValue[i].second);
            nValueRet += vValue[i].first;
        }
        return true;
    }

    if (nTotalLower < nTargetValue)
    {
        if (coinLowestLarger.second.first == NULL) {
            if (fLogFailureReason)
                LogPrintf("CWallet::SelectCoinsMinConf failed: sum of lower outputs (%s) < target (%s), and no other outputs\n",
                          FormatMoney(nTotalLower), FormatMoney(nTargetValue));
            return false;
        }
        setCoinsRet.insert(coinLowestLarger.second);
        nValueRet += coinLowestLarger.first;
        return true;
    }

    // Solve subset sum by stochastic approximation
    sort(vValue.rbegin(), vValue.rend(), CompareValueOnly());
    vector<char> vfBest;
    int64_t nBest;
    int nSize;

    ApproximateBestSubset(vValue, nTotalLower, nTargetValue, vfBest, nBest, nSize, 1000);
    if (nBest != nTargetValue && nTotalLower >= nTargetValue + CENT)
        ApproximateBestSubset(vValue, nTotalLower, nTargetValue + CENT, vfBest, nBest, nSize, 1000);

    // If we have a bigger coin and the stochastic approximation found more than 500 outputs to spend,
    // tidy up the dust in the wallet by spending just 500 of the found outputs plus the bigger coin
    if (coinLowestLarger.second.first && nSize > nMaxOutputsToSpend)
    {
        setCoinsRet.insert(coinLowestLarger.second);
        nValueRet += coinLowestLarger.first;

        nSize = nMaxOutputsToSpend;
        for (int i = vValue.size() - 1; i >= 0; i--)
            if (vfBest[i])
            {
                if (nSize-- > 0) {
                    setCoinsRet.insert(vValue[i].second);
                    nValueRet += vValue[i].first;
                } else
                    vfBest[i] = false;
            }

        LogPrint("selectcoins", "SelectCoins() smallest bigger plus truncated subset: ");
        LogPrint("selectcoins", "%s ", FormatMoney(coinLowestLarger.first));
        for (unsigned int i = 0; i < vValue.size(); i++)
            if (vfBest[i])
                LogPrint("selectcoins", "%s ", FormatMoney(vValue[i].first));
        LogPrint("selectcoins", "total %s\n", FormatMoney(nValueRet));
    // If we have a bigger coin and (either the stochastic approximation didn't find a good solution,
    //                                   or the next bigger coin is closer), return the bigger coin
    } else
    if (coinLowestLarger.second.first &&
        ((nBest != nTargetValue && nBest < nTargetValue + CENT) || coinLowestLarger.first <= nBest))
    {
        setCoinsRet.insert(coinLowestLarger.second);
        nValueRet += coinLowestLarger.first;
    }
    else {
        for (unsigned int i = 0; i < vValue.size(); i++)
            if (vfBest[i])
            {
                setCoinsRet.insert(vValue[i].second);
                nValueRet += vValue[i].first;
            }

        LogPrint("selectcoins", "SelectCoins() best subset: ");
        for (unsigned int i = 0; i < vValue.size(); i++)
            if (vfBest[i])
                LogPrint("selectcoins", "%s ", FormatMoney(vValue[i].first));
        LogPrint("selectcoins", "total %s\n", FormatMoney(nBest));
    }

    return true;
}

bool CWallet::SelectCoins(int64_t nTargetValue, unsigned int nSpendTime, set<pair<const CWalletTx*,unsigned int> >& setCoinsRet, int64_t& nValueRet, const CCoinControl* coinControl) const
{
    vector<COutput> vCoins;
    AvailableCoins(vCoins, true, coinControl);

    // coin control -> return all selected outputs (we want all selected to go into the transaction for sure)
    if (coinControl && coinControl->HasSelected())
    {
        BOOST_FOREACH(const COutput& out, vCoins)
        {
            nValueRet += out.tx->vout[out.i].nValue;
            setCoinsRet.insert(make_pair(out.tx, out.i));
        }

        if (nValueRet < nTargetValue) {
            LogPrintf("CWallet::SelectCoins failed: coin control selects less value than required (%s < %s)\n",
                      FormatMoney(nValueRet), FormatMoney(nTargetValue));
            return false;
        }

        return true;
    }

    boost::function<bool (const CWallet*, int64_t, unsigned int, int, int, std::vector<COutput>, std::set<std::pair<const CWalletTx*,unsigned int> >&, int64_t&, bool)> f = fMinimizeCoinAge ? &CWallet::SelectCoinsMinConfByCoinAge : &CWallet::SelectCoinsMinConf;

    if (setSpendLastAddresses.size()) {
        int64_t nMiscValue = 0;
        vector<COutput> vMiscCoins;
        BOOST_FOREACH(COutput& output, vCoins) {
            CTxDestination address;
            const CTxOut &txout = output.tx->vout[output.i];
            if (!ExtractDestination(txout.scriptPubKey, address) || !setSpendLastAddresses.count(address)) {
                nMiscValue += txout.nValue;
                vMiscCoins.push_back(output);
            }
        }

        LogPrintf("we have %d (of %d) non-spendlast outputs valued as %s\n", int(vMiscCoins.size()), int(vCoins.size()), FormatMoney(nMiscValue));

        if (nMiscValue >= nTargetValue &&
            (f(this, nTargetValue, nSpendTime, 1, 10, vMiscCoins, setCoinsRet, nValueRet, false) ||
             f(this, nTargetValue, nSpendTime, 1,  1, vMiscCoins, setCoinsRet, nValueRet, false) ||
             f(this, nTargetValue, nSpendTime, 0,  1, vMiscCoins, setCoinsRet, nValueRet, false))) {
            LogPrintf("we can make the transaction without spending any -spendlast coins\n");
            return true;
        }
    }

    if (f(this, nTargetValue, nSpendTime, 1, 10, vCoins, setCoinsRet, nValueRet, false) ||
        f(this, nTargetValue, nSpendTime, 1,  1, vCoins, setCoinsRet, nValueRet, false) ||
        f(this, nTargetValue, nSpendTime, 0,  1, vCoins, setCoinsRet, nValueRet, true)) // this is our last chance - if it fails log the reason why
        return true;

    LogPrintf("CWallet::SelectCoins failed: no way to select coins worth %s\n", FormatMoney(nTargetValue));

    return false;
}

// Select some coins without random shuffle or best subset approximation
bool CWallet::SelectCoinsForStaking(int64_t nTargetValue, unsigned int nSpendTime, set<pair<const CWalletTx*,unsigned int> >& setCoinsRet, int64_t& nValueRet) const
{
    vector<COutput> vCoins;
    AvailableCoinsForStaking(vCoins, nSpendTime);

    setCoinsRet.clear();
    nValueRet = 0;

    BOOST_FOREACH(COutput output, vCoins)
    {
        const CWalletTx *pcoin = output.tx;
        int i = output.i;

        // Stop if we've chosen enough inputs
        if (nValueRet >= nTargetValue)
            break;

        int64_t n = pcoin->vout[i].nValue;

        pair<int64_t,pair<const CWalletTx*,unsigned int> > coin = make_pair(n,make_pair(pcoin, i));

        if (n >= nTargetValue)
        {
            // If input value is greater or equal to target then simply insert
            //    it into the current subset and exit
            setCoinsRet.insert(coin.second);
            nValueRet += coin.first;
            break;
        }
        else if (n < nTargetValue + CENT)
        {
            setCoinsRet.insert(coin.second);
            nValueRet += coin.first;
        }
    }

    return true;
}

bool CWallet::CreateTransaction(const vector<pair<CScript, int64_t> >& vecSend, CWalletTx& wtxNew, CReserveKey& reservekey, int64_t& nFeeRet, std::string strCLAMSpeech, const CCoinControl* coinControl)
{
    int64_t nValue = 0;
    BOOST_FOREACH (const PAIRTYPE(CScript, int64_t)& s, vecSend)
    {
        if (nValue < 0) {
            LogPrintf("CWallet::CreateTransaction failed: negative value\n");
            return false;
        }
        nValue += s.second;
    }
    if (vecSend.empty() || nValue < 0) {
        LogPrintf("CWallet::CreateTransaction failed: no inputs or negative total\n");
        return false;
    }

    wtxNew.BindWallet(this);

    // set clamSpeech when creating a transaction
    if (strCLAMSpeech.empty() && !(mapArgs["-clamspeech"] == "off"))
        strCLAMSpeech = GetDefaultClamSpeech();

    wtxNew.strCLAMSpeech = strCLAMSpeech;
    if (wtxNew.strCLAMSpeech.length() > MAX_TX_COMMENT_LEN)
        wtxNew.strCLAMSpeech.resize(MAX_TX_COMMENT_LEN);

    {
        LOCK2(cs_main, cs_wallet);
        {
            nFeeRet = nTransactionFee;
            LogPrint("fee", "[FEE] start with fee = %s\n", FormatMoney(nFeeRet));
            while (true)
            {
                wtxNew.vin.clear();
                wtxNew.vout.clear();
                wtxNew.fFromMe = true;

                int64_t nTotalValue = nValue + nFeeRet;
                double dPriority = 0;
                // vouts to the payees
                  BOOST_FOREACH (const PAIRTYPE(CScript, int64)& s, vecSend)
                {
                    CTxOut txout(s.second, s.first);
                    if (txout.IsDust())
                    {
                        LogPrintf("Transaction amount too small");
                        return false;
                    }
                    wtxNew.vout.push_back(txout);
                }
                         // Choose coins to use
                set<pair<const CWalletTx*,unsigned int> > setCoins;
                int64_t nValueIn = 0;
                if (!SelectCoins(nTotalValue, wtxNew.nTime, setCoins, nValueIn, coinControl))
                {
                    LogPrintf("Insufficient funds");
                    return false;
                }
                BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
                {
                    int64_t nCredit = pcoin.first->vout[pcoin.second].nValue;
                    dPriority += (double)nCredit * pcoin.first->GetDepthInMainChain();
                }

                int64_t nChange = nValueIn - nValue - nFeeRet;
                // if sub-cent change is required, the fee must be raised to at least MIN_TX_FEE
                // or until nChange becomes zero
                // NOTE: this depends on the exact behaviour of GetMinFee
                if (nFeeRet < MIN_TX_FEE && nChange > 0 && nChange < CENT)
                {
                    int64_t nMoveToFee = min(nChange, MIN_TX_FEE - nFeeRet);
                    LogPrint("fee", "[FEE] change %s is positive and sub-cent so moving %s to fee\n", FormatMoney(nChange), FormatMoney(nMoveToFee));
                    nChange -= nMoveToFee;
                    nFeeRet += nMoveToFee;
                    LogPrint("fee", "[FEE] new change %s and fee %s\n", FormatMoney(nChange), FormatMoney(nFeeRet));
                }

                if (nChange > 0)
                {
                    // Fill a vout to ourself
                    // TODO: pass in scriptChange instead of reservekey so
                    // change transaction isn't always pay-to-bitcoin-address
                    CScript scriptChange;

                    // coin control: send change to custom address
                    if (coinControl && !boost::get<CNoDestination>(&coinControl->destChange))
                        scriptChange.SetDestination(coinControl->destChange);

                    // send change to one of the specified change addresses, if specified at init
                    else if (vChangeAddresses.size())
                    {
                        CKeyID keyID = vChangeAddresses[GetRandInt(vChangeAddresses.size())];
                        scriptChange.SetDestination(keyID);
                    }

                    // no coin control: send change to newly generated address
                    else
                    {
                        // Note: We use a new key here to keep it from being obvious which side is the change.
                        //  The drawback is that by not reusing a previous key, the change may be lost if a
                        //  backup is restored, if the backup doesn't have the new private key for the change.
                        //  If we reused the old key, it would be possible to add code to look for and
                        //  rediscover unknown transactions that were written with keys of ours to recover
                        //  post-backup change.

                        // Reserve a new key pair from key pool
                        CPubKey vchPubKey;
                        assert(reservekey.GetReservedKey(vchPubKey)); // should never fail, as we just unlocked

                        scriptChange.SetDestination(vchPubKey.GetID());
                    }

                    // Insert change txn at random position:
                    vector<CTxOut>::iterator position = wtxNew.vout.begin()+GetRandInt(wtxNew.vout.size()+1);
                    wtxNew.vout.insert(position, CTxOut(nChange, scriptChange));
                }
                else
                    reservekey.ReturnKey();

                // Fill vin
                BOOST_FOREACH(const PAIRTYPE(const CWalletTx*,unsigned int)& coin, setCoins)
                    wtxNew.vin.push_back(CTxIn(coin.first->GetHash(),coin.second));

                // Sign
                int nIn = 0;
                BOOST_FOREACH(const PAIRTYPE(const CWalletTx*,unsigned int)& coin, setCoins)
                    if (!SignSignature(*this, *coin.first, wtxNew, nIn++)) {
                        LogPrintf("CWallet::CreateTransaction failed: signing failed\n");
                        return false;
                    }

                // Limit size
                unsigned int nBytes = ::GetSerializeSize(*(CTransaction*)&wtxNew, SER_NETWORK, PROTOCOL_VERSION);
                if (nBytes >= MAX_STANDARD_TX_SIZE) {
                    LogPrintf("CWallet::CreateTransaction failed: transaction too big (%d >= %d)\n", nBytes, MAX_STANDARD_TX_SIZE);
                    return false;
                }
                dPriority /= nBytes;

                // Check that enough fee is included
                int64_t nPayFee = nTransactionFee * (1 + (int64_t)nBytes / 1000);
                int64_t nMinFee = wtxNew.GetMinFee(1, false, GMF_SEND);

                LogPrint("fee", "[FEE] tx is %d bytes; nPayFee = %s; nMinFee = %s; current fee = %s\n",
                         nBytes, FormatMoney(nPayFee), FormatMoney(nMinFee), FormatMoney(nFeeRet));
                if (nFeeRet < max(nPayFee, nMinFee))
                {
                    nFeeRet = max(nPayFee, nMinFee);
                    LogPrint("fee", "[FEE] so we increased fee to %s and loop again\n", FormatMoney(nFeeRet));
                    continue;
                }

                LogPrint("fee", "[FEE] so final fee is %s\n", FormatMoney(nFeeRet));

                // Fill vtxPrev by copying from previous transactions vtxPrev
                wtxNew.AddSupportingTransactions();
                wtxNew.fTimeReceivedIsTxTime = true;

                break;
            }
        }
    }
    return true;
}

bool CWallet::CreateTransaction(CScript scriptPubKey, int64_t nValue, int64_t nCount, CWalletTx& wtxNew, CReserveKey& reservekey, int64_t& nFeeRet, std::string strCLAMSpeech, const CCoinControl* coinControl)
{
    vector< pair<CScript, int64_t> > vecSend;
    while (nCount--)
        vecSend.push_back(make_pair(scriptPubKey, nValue));
    return CreateTransaction(vecSend, wtxNew, reservekey, nFeeRet, strCLAMSpeech, coinControl);
}

bool CWallet::GetStakeWeight(uint64_t& nWeight)
{
    // Choose coins to use
    int64_t nBalance = GetBalance();

    if (nBalance <= nReserveBalance)
        return false;

    vector<const CWalletTx*> vwtxPrev;

    set<pair<const CWalletTx*,unsigned int> > setCoins;
    int64_t nValueIn = 0;

    if (!SelectCoinsForStaking(nBalance - nReserveBalance, GetTime(), setCoins, nValueIn))
        return false;

    if (setCoins.empty())
        return false;

    nWeight = 0;

    int64_t nCurrentTime = GetTime();

    LOCK2(cs_main, cs_wallet);
    BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
    {
        CTransaction tx;
        uint256 hashBlock = 0;
        if (!GetTransaction(pcoin.first->GetHash(), tx, hashBlock, true))
            continue;

        if (nCurrentTime - pcoin.first->nTime > nStakeMinAge)
            nWeight += pcoin.first->vout[pcoin.second].nValue;
    } 

   return true;
}

bool CWallet::GetExpectedStakeTime(uint64_t& nExpected)
{
    unsigned int nBits = GetNextTargetRequired(pindexBest, true);
    CBigNum bnTarget;
    bnTarget.SetCompact(nBits);
    
    int64_t nBalance = GetBalance();
    if (nBalance <= nReserveBalance)
        return false;
    
    set<pair<const CWalletTx*,unsigned int> > setCoins;
    int64_t nValueIn = 0;
    double fail = 1;

    if (!SelectCoinsForStaking(nBalance - nReserveBalance, GetTime(), setCoins, nValueIn))
        return false;

    if (setCoins.empty())
        return false;

    BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
    {
        CTransaction tx;
        uint256 hashBlock = 0;
        {
            LOCK2(cs_main, cs_wallet);
            if (!GetTransaction(pcoin.first->GetHash(), tx, hashBlock, true))
                continue;
        }

        // p(A or B) = p(not((not A) and (not B))) = 1 - (p(1 - A) * p(1 - B))
        int64_t nValue = pcoin.first->vout[pcoin.second].nValue;
        if (nValue)
            fail *= 1 - 1.0/(CBigNum(~uint256(0)) / (bnTarget * nValue)).getulong();
    }

    nExpected = (STAKE_TIMESTAMP_MASK + 1) / (1 - fail);
    return true;
}

bool CWallet::CreateCoinStake(const CKeyStore& keystore, unsigned int nBits, int64_t nSearchInterval, int64_t nFees, CTransaction& txNew, CKey& key)
{
    CBlockIndex* pindexPrev = pindexBest;
    CBigNum bnTargetPerCoinDay;
    bnTargetPerCoinDay.SetCompact(nBits);

    txNew.vin.clear();
    txNew.vout.clear();

    // Mark coin stake transaction
    CScript scriptEmpty;
    scriptEmpty.clear();
    txNew.vout.push_back(CTxOut(0, scriptEmpty)); // this creates txNew.vout[0]

    // Choose coins to use
    int64_t nBalance = GetBalance();

    LogPrint("stake", "\n\n[STAKE] searching for staking opportunity:\n");

    if (nBalance <= nReserveBalance) {
        LogPrint("stake", "[STAKE] balance %s <= reserve %s\n", FormatMoney(nBalance), FormatMoney(nReserveBalance));
        return false;
    }

    vector<const CWalletTx*> vwtxPrev;

    set<pair<const CWalletTx*,unsigned int> > setCoins;
    int64_t nValueIn = 0;

    // Select coins with suitable depth
    if (!SelectCoinsForStaking(nBalance - nReserveBalance, txNew.nTime, setCoins, nValueIn)) {
        LogPrint("stake", "[STAKE] failed to select coins for staking\n");
        return false;
    }

    if (setCoins.empty()) {
        LogPrint("stake", "[STAKE] set of selected coins for staking is empty\n");
        return false;
    }

    LogPrint("stake", "[STAKE] checking %d output(s)\n", setCoins.size());

    int64_t nCredit = 0;
    CScript scriptPubKeyKernel;
    int64_t nBlockTime;
    CKeyID stakingkeyID;
    uint256 hashProofOfStake = 0;
    BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
    {
        if (!pcoin.first->hash)
            pcoin.first->hash = pcoin.first->GetHash();
        
        CTransaction tx;
        uint256 hashBlock = 0;
        {
            LOCK2(cs_main, cs_wallet);
            if (!GetTransaction(pcoin.first->hash, tx, hashBlock, true))
                continue;
        }

        // Read block header
        if (!mapBlockIndex.count(hashBlock))
            continue;

        CBlock block;
        {
            LOCK2(cs_main, cs_wallet);
            if (!block.ReadFromDisk(mapBlockIndex[hashBlock]))
                continue;
        }

            

        static int nMaxStakeSearchInterval = 60;
        nBlockTime = block.GetBlockTime();
        if (nBlockTime + nStakeMinAge > txNew.nTime - nMaxStakeSearchInterval) {
            LogPrint("stake", "[STAKE] skip %s:%-3d (%s CLAM) - only %d minutes old\n",
                      pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue),
                      (txNew.nTime - nMaxStakeSearchInterval - nBlockTime) / 60);
            continue; // only count coins meeting min age requirement
        }

        bool fKernelFound = false;
        for (unsigned int n=0; n<min(nSearchInterval,(int64_t)nMaxStakeSearchInterval) && !fKernelFound && pindexPrev == pindexBest; n++)
        {
            boost::this_thread::interruption_point();
            // Search backward in time from the given txNew timestamp 
            // Search nSearchInterval seconds back up to nMaxStakeSearchInterval
            uint256 targetProofOfStake = 0;
            COutPoint prevoutStake = COutPoint(pcoin.first->hash, pcoin.second);
            LogPrint("stake", "[STAKE] check %s:%-3d (%s CLAM)\n",
                      pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue));
            if (CheckStakeKernelHash(pindexPrev, nBits, block, prevoutStake.n, *pcoin.first, prevoutStake, txNew.nTime - n, hashProofOfStake, targetProofOfStake))
            {
                // Found a kernel
                LogPrint("coinstake", "CreateCoinStake : kernel found\n");
                vector<valtype> vSolutions;
                txnouttype whichType;
                CScript scriptPubKeyOut;
                scriptPubKeyKernel = pcoin.first->vout[pcoin.second].scriptPubKey;
                if (!Solver(scriptPubKeyKernel, whichType, vSolutions))
                {
                    LogPrint("stake", "[STAKE] fail %s:%-3d (%s CLAM) - can't parse kernel\n",
                              pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue));
                    break;
                }
                LogPrint("coinstake", "CreateCoinStake : parsed kernel type=%d\n", whichType);
                if (whichType != TX_PUBKEY && whichType != TX_PUBKEYHASH)
                {
                    LogPrint("stake", "[STAKE] fail %s:%-3d (%s CLAM) - bad kernel type\n",
                              pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue));
                    break;  // only support pay to public key and pay to address
                }
                if (whichType == TX_PUBKEYHASH) // pay to address type
                {
                    // convert to pay to public key type
                    stakingkeyID = uint160(vSolutions[0]);
                    if (!keystore.GetKey(stakingkeyID, key))
                    {
                        LogPrint("stake", "[STAKE] fail %s:%-3d (%s CLAM) - can't get public key (a)\n",
                                  pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue));
                        break;  // unable to find corresponding public key
                    }
                    scriptPubKeyOut << key.GetPubKey() << OP_CHECKSIG;
                }
                if (whichType == TX_PUBKEY)
                {
                    valtype& vchPubKey = vSolutions[0];
                    stakingkeyID = Hash160(vchPubKey);
                    if (!keystore.GetKey(stakingkeyID, key))
                    {
                        LogPrint("stake", "[STAKE] fail %s:%-3d (%s CLAM) - can't get public key, type %d\n",
                                  pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue),
                                  whichType);
                        LogPrint("coinstake", "CreateCoinStake : failed to get key for kernel type=%d\n", whichType);
                        break;  // unable to find corresponding public key
                    }

                    if (key.GetPubKey() != vchPubKey)
                    {
                        LogPrint("stake", "[STAKE] fail %s:%-3d (%s CLAM) - invalid key, type %d\n",
                                  pcoin.first->hash.ToString(), pcoin.second, FormatMoney(pcoin.first->vout[pcoin.second].nValue),
                                  whichType);
                        break; // keys mismatch
                    }

                    scriptPubKeyOut = scriptPubKeyKernel;
                }

                txNew.nTime -= n;
                txNew.vin.push_back(CTxIn(pcoin.first->hash, pcoin.second));
                nCredit += pcoin.first->vout[pcoin.second].nValue;
                vwtxPrev.push_back(pcoin.first);
                txNew.vout.push_back(CTxOut(0, scriptPubKeyOut)); // this creates txNew.vout[1]
                setCoins.erase(pcoin); // don't consider the staking coin for merging later

                LogPrint("coinstake", "CreateCoinStake : added kernel type=%d\n", whichType);
                fKernelFound = true;
                break;
            }
        }

        if (fKernelFound)
            break; // if kernel is found stop searching
    }

    if (nCredit == 0 || nCredit > nBalance - nReserveBalance)
        return false;

    // Calculate coin age reward
    int64_t nReward;
    {
        uint64_t nCoinAge;
        if (!txNew.GetCoinAge(nCoinAge))
            return error("CreateCoinStake : failed to calculate coin age");

        nReward = GetProofOfStakeReward(pindexBest, nCoinAge, nFees);
        if (nReward <= 0)
            return false;

        if (!fRewardTo)
            nCredit += nReward;
    }

    // set clamSpeech when staking a block
    if (!(mapArgs["-clamstake"] == "off")) {
        if (weightedStakeSpeech.size()) {
            txNew.strCLAMSpeech = weightedStakeSpeech.select(hashProofOfStake.Get64());
        } else if(GetDefaultClamourClamSpeech() == "") {
            txNew.strCLAMSpeech = GetDefaultClamSpeech();
        } else {
            txNew.strCLAMSpeech = GetDefaultClamourClamSpeech();
        }
        if (txNew.strCLAMSpeech.length() > MAX_TX_COMMENT_LEN)
            txNew.strCLAMSpeech.resize(MAX_TX_COMMENT_LEN);
    }

    CScript scriptStakeTo;
    int nFirstRealOutput;
    // if we are staking to a specific address, and it's different than the address that is staking, make a new output for it
    if (fStakeTo && stakingkeyID != staketokeyID) {
        nFirstRealOutput = 2;
        scriptStakeTo.SetDestination(staketokeyID);
        txNew.vout.push_back(CTxOut(0, scriptStakeTo));
    } else {
        nFirstRealOutput = 1;
        scriptStakeTo = txNew.vout[1].scriptPubKey;
    }

    // if we're splitting on size, only split if the new size will make it at least double the split size
    // and if we're not, split on age
    if (( nSplitSize && nCredit >= nSplitSize * 2) ||
        (!nSplitSize && GetWeight(nBlockTime, (int64_t)txNew.nTime) < GetStakeSplitAge()))
    {
        if (nSplitSize) {
            int n = nFirstRealOutput;

            while (true) {
                if (n > nFirstRealOutput)
                    txNew.vout.push_back(CTxOut(0, scriptStakeTo));

                if (nCredit < nSplitSize * 2) {
                    txNew.vout[n].nValue = nCredit;
                    break;
                }

                txNew.vout[n].nValue = nSplitSize;
                nCredit -= nSplitSize;

                n++;
            }
        } else {
            txNew.vout.push_back(CTxOut(0, scriptStakeTo)); //split stake
            txNew.vout[nFirstRealOutput].nValue = (nCredit / 2 / CENT) * CENT;
            txNew.vout[nFirstRealOutput+1].nValue = nCredit - txNew.vout[nFirstRealOutput].nValue;
        }
    } else {
        // we're not splitting the output, so attempt to add more inputs
        BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
        {
            // Only add coins of the same key/address as kernel, or if fCombineAny is true
            if ((fCombineAny ||
                 pcoin.first->vout[pcoin.second].scriptPubKey == scriptPubKeyKernel ||
                 pcoin.first->vout[pcoin.second].scriptPubKey == txNew.vout[1].scriptPubKey) &&
                (pcoin.first->hash != txNew.vin[0].prevout.hash ||
                 pcoin.second != txNew.vin[0].prevout.n))
            {
                // Only add coins of the same key/address as kernel, or if fCombineAny is true
                if (fCombineAny ||
                    pcoin.first->vout[pcoin.second].scriptPubKey == scriptPubKeyKernel ||
                    pcoin.first->vout[pcoin.second].scriptPubKey == txNew.vout[1].scriptPubKey)
                {
                    // skip this output if ..
                    if (nCredit + pcoin.first->vout[pcoin.second].nValue > nCombineLimit) // .. it causes the total to exceed the combine threshold
                        continue;

                    // skip all remaining outputs if..
                    if (txNew.vin.size() >= 100 ||                                                      // .. we already have enough inputs
                        nCredit + pcoin.first->vout[pcoin.second].nValue > nBalance - nReserveBalance)  // or we have reached the reserve limit
                        break;

                txNew.vin.push_back(CTxIn(pcoin.first->hash, pcoin.second));
                nCredit += pcoin.first->vout[pcoin.second].nValue;
                vwtxPrev.push_back(pcoin.first);
            }
        }

        // we're not splitting, so put all the value in the first and only real output
        txNew.vout[nFirstRealOutput].nValue = nCredit;
    }

    // if we're sending just the staking reward to a separate address, add that output on at the end
    if (fRewardTo) {
        scriptStakeTo.SetDestination(rewardtokeyID);
        txNew.vout.push_back(CTxOut(nReward, scriptStakeTo));
    }

    // Sign
    int nIn = 0;
    BOOST_FOREACH(const CWalletTx* pcoin, vwtxPrev)
    {
        if (!SignSignature(*this, *pcoin, txNew, nIn++))
            return error("CreateCoinStake : failed to sign coinstake");
    }

    // Limit size
    unsigned int nBytes = ::GetSerializeSize(txNew, SER_NETWORK, PROTOCOL_VERSION);
    if (nBytes >= MAX_BLOCK_SIZE_GEN/5)
        return error("CreateCoinStake : exceeded coinstake size limit");

    // Successfully generated coinstake
    return true;
}


// Call after CreateTransaction unless you want to abort
bool CWallet::CommitTransaction(CWalletTx& wtxNew, CReserveKey& reservekey)
{
    {
        LOCK2(cs_main, cs_wallet);
        LogPrintf("CommitTransaction:\n%s", wtxNew.ToString());
        {
            // This is only to keep the database open to defeat the auto-flush for the
            // duration of this scope.  This is the only place where this optimization
            // maybe makes sense; please don't do it anywhere else.
            CWalletDB* pwalletdb = fFileBacked ? new CWalletDB(strWalletFile,"r") : NULL;

            // Take key pair from key pool so it won't be used again
            reservekey.KeepKey();

            // Add tx to wallet, because if it has change it's also ours,
            // otherwise just for transaction history.
            AddToWallet(wtxNew);

            // Mark old coins as spent
            set<CWalletTx*> setCoins;
            BOOST_FOREACH(const CTxIn& txin, wtxNew.vin)
            {
                CWalletTx &coin = mapWallet[txin.prevout.hash];
                coin.BindWallet(this);
                coin.MarkSpent(txin.prevout.n);
                coin.WriteToDisk();
                NotifyTransactionChanged(this, coin.GetHash(), CT_UPDATED);
            }

            if (fFileBacked)
                delete pwalletdb;
        }

        // Track how many getdata requests our transaction gets
        mapRequestCount[wtxNew.GetHash()] = 0;

        // Broadcast
        if (!wtxNew.AcceptToMemoryPool(true))
        {
            // This must not fail. The transaction has already been signed and recorded.
            LogPrintf("CommitTransaction() : Error: Transaction not valid\n");
            return false;
        }
        wtxNew.RelayWalletTransaction();
    }
    return true;
}




string CWallet::SendMoney(CScript scriptPubKey, int64_t nValue, int64_t nCount, CWalletTx& wtxNew, std::string strTxComment, bool fAskFee)
{
    CReserveKey reservekey(this);
    int64_t nFeeRequired;

    if (IsLocked())
    {
        string strError = _("Error: Wallet locked, unable to create transaction!");
        LogPrintf("SendMoney() : %s", strError);
        return strError;
    }
    if (fWalletUnlockStakingOnly)
    {
        string strError = _("Error: Wallet unlocked for staking only, unable to create transaction.");
        LogPrintf("SendMoney() : %s", strError);
        return strError;
    }
    if (!CreateTransaction(scriptPubKey, nValue, nCount, wtxNew, reservekey, nFeeRequired, strTxComment))
    {
        string strError;
        if (nValue * nCount + nFeeRequired > GetBalance())
            strError = strprintf(_("Error: This transaction requires a transaction fee of at least %s because of its amount, complexity, or use of recently received funds!"), FormatMoney(nFeeRequired));
        else
            strError = _("Error: Transaction creation failed!");
        LogPrintf("SendMoney() : %s\n", strError);
        return strError;
    }

    if (fAskFee && !uiInterface.ThreadSafeAskFee(nFeeRequired, _("Sending...")))
        return "ABORTED";

    if (!CommitTransaction(wtxNew, reservekey))
        return _("Error: The transaction was rejected! This might happen if some of the coins in your wallet were already spent, such as if you used a copy of wallet.dat and coins were spent in the copy but not marked as spent here.");

    return "";
}

string CWallet::SendCLAMSpeech(CWalletTx& wtxNew, string clamSpeech, string prefix, bool fAskFee)
{ 
    if (prefix == "notary") 
    {
        uint256 hash;
        hash.SetHex(clamSpeech);
        clamSpeech = "notary " + hash.GetHex();

    } 
    else if (prefix == "clamour") 
    {
        clamSpeech = "create clamour" + clamSpeech;
    } 
    else if (prefix.length() > 0)
    {
        clamSpeech = prefix + clamSpeech;
    } 

    CReserveKey reservekey(this);
    int64_t nFeeRequired;

    if (GetBalance() <= 0 ) {
        LogPrintf("CWallet::SendNotary failed: you require a balance to post a notary entry\n",
                  FormatMoney(nTransactionFee), FormatMoney(GetBalance()));
        return _("Insufficient funds");
    }

    if (IsLocked())
    {
        string strError = _("Error: Wallet locked, unable to create norary transaction!");
        LogPrintf("SendNotary() : %s", strError);
        return strError;
    }
    if (fWalletUnlockStakingOnly)
    {
        string strError = _("Error: Wallet unlocked for staking only, unable to create notary transaction.");
        LogPrintf("SendNotary() : %s", strError);
        return strError;
    }
    if (!CreateCLAMSpeechTransaction(wtxNew, reservekey, nFeeRequired, clamSpeech))
    {
        string strError;
        if (nFeeRequired > GetBalance())
            strError = strprintf(_("Error: This notary transaction requires a transaction fee of at least %s because of its amount, complexity, or use of recently received funds!"), FormatMoney(nFeeRequired));
        else
            strError = _("Error: Notary transaction creation failed!");
        LogPrintf("SendNotary() : %s\n", strError);

        return strError;
    }

    if (fAskFee && !uiInterface.ThreadSafeAskFee(nFeeRequired, _("Sending...")))
        return "ABORTED";

    if (!CommitTransaction(wtxNew, reservekey))
        return _("Error: The transaction was rejected! This might happen if some of the coins in your wallet were already spent, such as if you used a copy of wallet.dat and coins were spent in the copy but not marked as spent here.");

    return "";


}

bool CWallet::CreateCLAMSpeechTransaction(CWalletTx& wtxNew, CReserveKey& reservekey, int64_t& nFeeRet, std::string clamSpeech, const CCoinControl* coinControl) 
{
	//vector< pair<CScript, int64_t> > vecSend;
	//vecSend.push_back(make_pair(scriptPubKey, nValue));

	wtxNew.BindWallet(this);

    wtxNew.strCLAMSpeech = clamSpeech;   

	{
		LOCK2(cs_main, cs_wallet);
		{
			nFeeRet = nTransactionFee;
			LogPrint("fee", "[FEE] start with fee = %s\n", FormatMoney(nFeeRet));
			while(true)
			{
				wtxNew.vin.clear();
				wtxNew.vout.clear();
				wtxNew.fFromMe = true;

				double dPriority = 0;
				
				set<pair<const CWalletTx*,unsigned int> > setCoins;
				int64_t nValueIn=0;
				if(!SelectCoins(1, wtxNew.nTime, setCoins, nValueIn, coinControl)) {
					LogPrintf("CWallet::CreateNotaryTransaction failed: coin selection failed\n");
				}

				BOOST_FOREACH(PAIRTYPE(const CWalletTx*, unsigned int) pcoin, setCoins)
                		{
                    			int64_t nCredit = pcoin.first->vout[pcoin.second].nValue;
                    			dPriority += (double)nCredit * pcoin.first->GetDepthInMainChain();
               		 	}

				int64_t nChange = nValueIn - nFeeRet;

				// NOTE: this depends on the exact behaviour of GetMinFee
                		if (nFeeRet < MIN_TX_FEE && nChange > 0 && nChange < CENT)
                		{
                    			int64_t nMoveToFee = min(nChange, MIN_TX_FEE - nFeeRet);
                    			LogPrint("fee", "[FEE] change %s is positive and sub-cent so moving %s to fee\n", FormatMoney(nChange), FormatMoney(nMoveToFee));
                    			nChange -= nMoveToFee;
                    			nFeeRet += nMoveToFee;
                    			LogPrint("fee", "[FEE] new change %s and fee %s\n", FormatMoney(nChange), FormatMoney(nFeeRet));
                		}
				
				CScript scriptNotary; 
				CPubKey vchPubKey;
				assert(reservekey.GetReservedKey(vchPubKey));

				scriptNotary.SetDestination(vchPubKey.GetID());

               			wtxNew.vout.push_back(CTxOut(nChange, scriptNotary));

		                // Fill vin
                		BOOST_FOREACH(const PAIRTYPE(const CWalletTx*,unsigned int)& coin, setCoins)
                   		wtxNew.vin.push_back(CTxIn(coin.first->GetHash(),coin.second));

				// Sign
               	 		int nIn = 0;
                		BOOST_FOREACH(const PAIRTYPE(const CWalletTx*,unsigned int)& coin, setCoins)
                    		if (!SignSignature(*this, *coin.first, wtxNew, nIn++)) {
                        		LogPrintf("CWallet::CreateTransaction failed: signing failed\n");
                        		return false;
                    		}

                		// Limit size
                		unsigned int nBytes = ::GetSerializeSize(*(CTransaction*)&wtxNew, SER_NETWORK, PROTOCOL_VERSION);
                		if (nBytes >= MAX_STANDARD_TX_SIZE) {
                    			LogPrintf("CWallet::CreateTransaction failed: transaction too big (%d >= %d)\n", nBytes, MAX_STANDARD_TX_SIZE);
                    			return false;
                		}
                		dPriority /= nBytes;

                		// Check that enough fee is included
                		int64_t nPayFee = nTransactionFee * (1 + (int64_t)nBytes / 1000);
                		int64_t nMinFee = wtxNew.GetMinFee(1, false, GMF_SEND);

                			LogPrint("fee", "[FEE] tx is %d bytes; nPayFee = %s; nMinFee = %s; current fee = %s\n",
                        		nBytes, FormatMoney(nPayFee), FormatMoney(nMinFee), FormatMoney(nFeeRet));
                		if (nFeeRet < max(nPayFee, nMinFee))
                		{
                   	 		nFeeRet = max(nPayFee, nMinFee);
                    			LogPrint("fee", "[FEE] so we increased fee to %s and loop again\n", FormatMoney(nFeeRet));
                    			continue;
                		}

                		LogPrint("fee", "[FEE] so final fee is %s\n", FormatMoney(nFeeRet));

                		// Fill vtxPrev by copying from previous transactions vtxPrev
                		wtxNew.AddSupportingTransactions();
                		wtxNew.fTimeReceivedIsTxTime = true;

                		break;
			}
		}
	}
	return true;
}

string CWallet::SendMoneyToDestination(const CTxDestination& address, int64_t nValue, int64_t nCount, CWalletTx& wtxNew, std::string strTxComment, bool fAskFee)
{
    // Check amount
    if (nValue <= 0)
        return _("Invalid amount");
    // Check count
    if (nCount <= 0)
        return _("Invalid count");
    int64_t nValueProduct = nValue * nCount;
    if (nValueProduct + nTransactionFee > GetBalance()) {
        LogPrintf("CWallet::SendMoneyToDestination failed: value (%s) + fee (%s) = %s > balance (%s)\n",
                  FormatMoney(nValueProduct), FormatMoney(nTransactionFee), FormatMoney(nValueProduct + nTransactionFee), FormatMoney(GetBalance()));
        return _("Insufficient funds");
    }

    // Parse Bitcoin address
    CScript scriptPubKey;
    scriptPubKey.SetDestination(address);

    return SendMoney(scriptPubKey, nValue, nCount, wtxNew, strTxComment, fAskFee);
}




DBErrors CWallet::LoadWallet(bool& fFirstRunRet)
{
    if (!fFileBacked)
        return DB_LOAD_OK;
    fFirstRunRet = false;
    DBErrors nLoadWalletRet = CWalletDB(strWalletFile,"cr+").LoadWallet(this);
    if (nLoadWalletRet == DB_NEED_REWRITE)
    {
        if (CDB::Rewrite(strWalletFile, "\x04pool"))
        {
            LOCK(cs_wallet);
            setKeyPool.clear();
            // Note: can't top-up keypool here, because wallet is locked.
            // User will be prompted to unlock wallet the next operation
            // the requires a new key.
        }
    }

    if (nLoadWalletRet != DB_LOAD_OK)
        return nLoadWalletRet;
    fFirstRunRet = !vchDefaultKey.IsValid();

    return DB_LOAD_OK;
}

DBErrors CWallet::LoadWalletImport()
{
    if (!fFileBacked)
        return DB_LOAD_OK;
    DBErrors nLoadWalletRet = CWalletDB(strWalletFile,"r+").LoadWalletImport(this);
    if (nLoadWalletRet != DB_LOAD_OK) {
        return nLoadWalletRet;
    }
    
    return DB_LOAD_OK;
}


bool CWallet::SetAddressBookName(const CTxDestination& address, const string& strName)
{
    bool fUpdated = false;
    {
        LOCK(cs_wallet); // mapAddressBook
        std::map<CTxDestination, std::string>::iterator mi = mapAddressBook.find(address);
        fUpdated = mi != mapAddressBook.end();
        mapAddressBook[address] = strName;
    }
    NotifyAddressBookChanged(this, address, strName, ::IsMine(*this, address),
                             (fUpdated ? CT_UPDATED : CT_NEW) );
    if (!fFileBacked)
        return false;
    return CWalletDB(strWalletFile).WriteName(CBitcoinAddress(address).ToString(), strName);
}

bool CWallet::DelAddressBookName(const CTxDestination& address)
{
    {
        LOCK(cs_wallet); // mapAddressBook

        mapAddressBook.erase(address);
    }

    NotifyAddressBookChanged(this, address, "", ::IsMine(*this, address), CT_DELETED);

    if (!fFileBacked)
        return false;
    return CWalletDB(strWalletFile).EraseName(CBitcoinAddress(address).ToString());
}

bool CWallet::SetDefaultKey(const CPubKey &vchPubKey)
{
    if (fFileBacked)
    {
        if (!CWalletDB(strWalletFile).WriteDefaultKey(vchPubKey))
            return false;
    }
    vchDefaultKey = vchPubKey;
    return true;
}

//
// Mark old keypool keys as used,
// and generate all new keys
//
bool CWallet::NewKeyPool()
{
    {
        LOCK(cs_wallet);
        CWalletDB walletdb(strWalletFile);
        BOOST_FOREACH(int64_t nIndex, setKeyPool)
            walletdb.ErasePool(nIndex);
        setKeyPool.clear();

        if (IsLocked())
            return false;

        int64_t nKeys = max(GetArg("-keypool", 100), (int64_t)0);
        for (int i = 0; i < nKeys; i++)
        {
            int64_t nIndex = i+1;
            walletdb.WritePool(nIndex, CKeyPool(GenerateNewKey()));
            setKeyPool.insert(nIndex);
        }
        LogPrintf("CWallet::NewKeyPool wrote %d new keys\n", nKeys);
    }
    return true;
}

bool CWallet::TopUpKeyPool(unsigned int nSize)
{
    {
        LOCK(cs_wallet);

        if (IsLocked())
            return false;

        CWalletDB walletdb(strWalletFile);

        // Top up key pool
        unsigned int nTargetSize;
        if (nSize > 0)
            nTargetSize = nSize;
        else
            nTargetSize = max(GetArg("-keypool", 100), (int64_t)0);

        while (setKeyPool.size() < (nTargetSize + 1))
        {
            int64_t nEnd = 1;
            if (!setKeyPool.empty())
                nEnd = *(--setKeyPool.end()) + 1;
            if (!walletdb.WritePool(nEnd, CKeyPool(GenerateNewKey())))
                throw runtime_error("TopUpKeyPool() : writing generated key failed");
            setKeyPool.insert(nEnd);
            LogPrintf("keypool added key %d, size=%"PRIszu"\n", nEnd, setKeyPool.size());
        }
    }
    return true;
}

void CWallet::ReserveKeyFromKeyPool(int64_t& nIndex, CKeyPool& keypool)
{
    nIndex = -1;
    keypool.vchPubKey = CPubKey();
    {
        LOCK(cs_wallet);

        if (!IsLocked())
            TopUpKeyPool();

        // Get the oldest key
        if(setKeyPool.empty())
            return;

        CWalletDB walletdb(strWalletFile);

        nIndex = *(setKeyPool.begin());
        setKeyPool.erase(setKeyPool.begin());
        if (!walletdb.ReadPool(nIndex, keypool))
            throw runtime_error("ReserveKeyFromKeyPool() : read failed");
        if (!HaveKey(keypool.vchPubKey.GetID()))
            throw runtime_error("ReserveKeyFromKeyPool() : unknown key in key pool");
        assert(keypool.vchPubKey.IsValid());
        LogPrintf("keypool reserve %d\n", nIndex);
    }
}

int64_t CWallet::AddReserveKey(const CKeyPool& keypool)
{
    {
        LOCK2(cs_main, cs_wallet);
        CWalletDB walletdb(strWalletFile);

        int64_t nIndex = 1 + *(--setKeyPool.end());
        if (!walletdb.WritePool(nIndex, keypool))
            throw runtime_error("AddReserveKey() : writing added key failed");
        setKeyPool.insert(nIndex);
        return nIndex;
    }
    return -1;
}

void CWallet::KeepKey(int64_t nIndex)
{
    // Remove from key pool
    if (fFileBacked)
    {
        CWalletDB walletdb(strWalletFile);
        walletdb.ErasePool(nIndex);
    }
    LogPrintf("keypool keep %d\n", nIndex);
}

void CWallet::ReturnKey(int64_t nIndex)
{
    // Return to key pool
    {
        LOCK(cs_wallet);
        setKeyPool.insert(nIndex);
    }
    LogPrintf("keypool return %d\n", nIndex);
}

bool CWallet::GetKeyFromPool(CPubKey& result)
{
    int64_t nIndex = 0;
    CKeyPool keypool;
    {
        LOCK(cs_wallet);
        ReserveKeyFromKeyPool(nIndex, keypool);
        if (nIndex == -1)
        {
            if (IsLocked()) return false;
            result = GenerateNewKey();
            return true;
        }
        KeepKey(nIndex);
        result = keypool.vchPubKey;
    }
    return true;
}

int64_t CWallet::GetOldestKeyPoolTime()
{
    int64_t nIndex = 0;
    CKeyPool keypool;
    ReserveKeyFromKeyPool(nIndex, keypool);
    if (nIndex == -1)
        return GetTime();
    ReturnKey(nIndex);
    return keypool.nTime;
}

std::map<CTxDestination, int64_t> CWallet::GetAddressBalances()
{
    map<CTxDestination, int64_t> balances;

    {
        LOCK(cs_wallet);
        BOOST_FOREACH(PAIRTYPE(uint256, CWalletTx) walletEntry, mapWallet)
        {
            CWalletTx *pcoin = &walletEntry.second;

            if (!IsFinalTx(*pcoin) || !pcoin->IsTrusted())
                continue;

            if ((pcoin->IsCoinBase() || pcoin->IsCoinStake()) && pcoin->GetBlocksToMaturity() > 0)
                continue;

            int nDepth = pcoin->GetDepthInMainChain();
            if (nDepth < (pcoin->IsFromMe() ? 0 : 1))
                continue;

            for (unsigned int i = 0; i < pcoin->vout.size(); i++)
            {
                CTxDestination addr;
                if (!IsMine(pcoin->vout[i]))
                    continue;
                if(!ExtractDestination(pcoin->vout[i].scriptPubKey, addr))
                    continue;

                int64_t n = pcoin->IsSpent(i) ? 0 : pcoin->vout[i].nValue;

                if (!balances.count(addr))
                    balances[addr] = 0;
                balances[addr] += n;
            }
        }
    }

    return balances;
}

set< set<CTxDestination> > CWallet::GetAddressGroupings()
{
    AssertLockHeld(cs_wallet); // mapWallet
    set< set<CTxDestination> > groupings;
    set<CTxDestination> grouping;

    BOOST_FOREACH(PAIRTYPE(uint256, CWalletTx) walletEntry, mapWallet)
    {
        CWalletTx *pcoin = &walletEntry.second;

        if (pcoin->vin.size() > 0 && IsMine(pcoin->vin[0]))
        {
            // group all input addresses with each other
            BOOST_FOREACH(CTxIn txin, pcoin->vin)
            {
                CTxDestination address;
                if(!ExtractDestination(mapWallet[txin.prevout.hash].vout[txin.prevout.n].scriptPubKey, address))
                    continue;
                grouping.insert(address);
            }

            // group change with input addresses
            BOOST_FOREACH(CTxOut txout, pcoin->vout)
                if (IsChange(txout))
                {
                    CWalletTx tx = mapWallet[pcoin->vin[0].prevout.hash];
                    CTxDestination txoutAddr;
                    if(!ExtractDestination(txout.scriptPubKey, txoutAddr))
                        continue;
                    grouping.insert(txoutAddr);
                }
            groupings.insert(grouping);
            grouping.clear();
        }

        // group lone addrs by themselves
        for (unsigned int i = 0; i < pcoin->vout.size(); i++)
            if (IsMine(pcoin->vout[i]))
            {
                CTxDestination address;
                if(!ExtractDestination(pcoin->vout[i].scriptPubKey, address))
                    continue;
                grouping.insert(address);
                groupings.insert(grouping);
                grouping.clear();
            }
    }

    set< set<CTxDestination>* > uniqueGroupings; // a set of pointers to groups of addresses
    map< CTxDestination, set<CTxDestination>* > setmap;  // map addresses to the unique group containing it
    BOOST_FOREACH(set<CTxDestination> grouping, groupings)
    {
        // make a set of all the groups hit by this new group
        set< set<CTxDestination>* > hits;
        map< CTxDestination, set<CTxDestination>* >::iterator it;
        BOOST_FOREACH(CTxDestination address, grouping)
            if ((it = setmap.find(address)) != setmap.end())
                hits.insert((*it).second);

        // merge all hit groups into a new single group and delete old groups
        set<CTxDestination>* merged = new set<CTxDestination>(grouping);
        BOOST_FOREACH(set<CTxDestination>* hit, hits)
        {
            merged->insert(hit->begin(), hit->end());
            uniqueGroupings.erase(hit);
            delete hit;
        }
        uniqueGroupings.insert(merged);

        // update setmap
        BOOST_FOREACH(CTxDestination element, *merged)
            setmap[element] = merged;
    }

    set< set<CTxDestination> > ret;
    BOOST_FOREACH(set<CTxDestination>* uniqueGrouping, uniqueGroupings)
    {
        ret.insert(*uniqueGrouping);
        delete uniqueGrouping;
    }

    return ret;
}

// ppcoin: disable transaction (only for coinstake)
void CWallet::DisableTransaction(const CTransaction &tx)
{
    if (!tx.IsCoinStake() || !IsFromMe(tx))
        return; // only disconnecting coinstake requires marking input unspent

    LOCK(cs_wallet);
    BOOST_FOREACH(const CTxIn& txin, tx.vin)
    {
        map<uint256, CWalletTx>::iterator mi = mapWallet.find(txin.prevout.hash);
        if (mi != mapWallet.end())
        {
            CWalletTx& prev = (*mi).second;
            if (txin.prevout.n < prev.vout.size() && IsMine(prev.vout[txin.prevout.n]))
            {
                prev.MarkUnspent(txin.prevout.n);
                prev.WriteToDisk();
            }
        }
    }
}

void CWallet::SumStakingRewards()
{
    mapAddressRewards.clear();

    int nWalletTx = 0, nStakeTx = 0;

    string strAccount;
    int64_t nFee;
    CTxDestination address;
    list<pair<CTxDestination, int64_t> > listReceived, listSent;

    typedef std::map<uint256, CWalletTx> TMapWallet;
    BOOST_FOREACH(TMapWallet::value_type& entry, mapWallet)
    {
        const CWalletTx& wtx = entry.second;
        if (wtx.IsCoinStake() &&
            IsFinalTx(wtx) &&
            wtx.GetDepthInMainChain() > 0 &&
            wtx.vout.size() > 1) {

            if (!ExtractDestination(wtx.vout[1].scriptPubKey, address)) {
                LogPrintf("can't extract destination\n");
                continue;
            }

            std::string addr(CBitcoinAddress(address).ToString());
            // this can happen when the stake transaction has multiple outputs, and one of them goes to us
            if (!::IsMine(*this, address)) {
                // LogPrintf("staked address %s is not mine in txid %s\n", addr, wtx.GetHash().ToString());
                continue;
            }

            wtx.GetAmounts(listReceived, listSent, nFee, strAccount);
            mapAddressRewards["*"] -= nFee;
            mapAddressRewards[addr] -= nFee;
            nStakeTx++;
        }
        
        nWalletTx++;
    }
    LogPrintf("CWallet::SumStakingRewards() found %d of %d tx to be stakes\n", nStakeTx, nWalletTx);
    fAddressRewardsReady = true;
}

bool CReserveKey::GetReservedKey(CPubKey& pubkey)
{
    if (nIndex == -1)
    {
        CKeyPool keypool;
        pwallet->ReserveKeyFromKeyPool(nIndex, keypool);
        if (nIndex != -1)
            vchPubKey = keypool.vchPubKey;
        else {
            if (pwallet->vchDefaultKey.IsValid()) {
                LogPrintf("CReserveKey::GetReservedKey(): Warning: Using default key instead of a new key, top up your keypool!");
                vchPubKey = pwallet->vchDefaultKey;
            } else
                return false;
        }
    }
    assert(vchPubKey.IsValid());
    pubkey = vchPubKey;
    return true;
}

void CReserveKey::KeepKey()
{
    if (nIndex != -1)
        pwallet->KeepKey(nIndex);
    nIndex = -1;
    vchPubKey = CPubKey();
}

void CReserveKey::ReturnKey()
{
    if (nIndex != -1)
        pwallet->ReturnKey(nIndex);
    nIndex = -1;
    vchPubKey = CPubKey();
}

void CWallet::GetAllReserveKeys(set<CKeyID>& setAddress) const
{
    setAddress.clear();

    CWalletDB walletdb(strWalletFile);

    LOCK2(cs_main, cs_wallet);
    BOOST_FOREACH(const int64_t& id, setKeyPool)
    {
        CKeyPool keypool;
        if (!walletdb.ReadPool(id, keypool))
            throw runtime_error("GetAllReserveKeyHashes() : read failed");
        assert(keypool.vchPubKey.IsValid());
        CKeyID keyID = keypool.vchPubKey.GetID();
        if (!HaveKey(keyID))
            throw runtime_error("GetAllReserveKeyHashes() : unknown key in key pool");
        setAddress.insert(keyID);
    }
}

void CWallet::UpdatedTransaction(const uint256 &hashTx)
{
    {
        LOCK(cs_wallet);
        // Only notify UI if this transaction is in this wallet
        map<uint256, CWalletTx>::const_iterator mi = mapWallet.find(hashTx);
        if (mi != mapWallet.end())
            NotifyTransactionChanged(this, hashTx, CT_UPDATED);
    }
}

void CWallet::GetKeyBirthTimes(std::map<CKeyID, int64_t> &mapKeyBirth) const {
    AssertLockHeld(cs_wallet); // mapKeyMetadata
    mapKeyBirth.clear();

    // get birth times for keys with metadata
    for (std::map<CKeyID, CKeyMetadata>::const_iterator it = mapKeyMetadata.begin(); it != mapKeyMetadata.end(); it++)
        if (it->second.nCreateTime)
            mapKeyBirth[it->first] = it->second.nCreateTime;

    // map in which we'll infer heights of other keys
    CBlockIndex *pindexMax = FindBlockByHeight(std::max(0, nBestHeight - 144)); // the tip can be reorganised; use a 144-block safety margin
    std::map<CKeyID, CBlockIndex*> mapKeyFirstBlock;
    std::set<CKeyID> setKeys;
    GetKeys(setKeys);
    BOOST_FOREACH(const CKeyID &keyid, setKeys) {
        if (mapKeyBirth.count(keyid) == 0)
            mapKeyFirstBlock[keyid] = pindexMax;
    }
    setKeys.clear();

    // if there are no such keys, we're done
    if (mapKeyFirstBlock.empty())
        return;

    // find first block that affects those keys, if there are any left
    std::vector<CKeyID> vAffected;
    for (std::map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); it++) {
        // iterate over all wallet transactions...
        const CWalletTx &wtx = (*it).second;
        std::map<uint256, CBlockIndex*>::const_iterator blit = mapBlockIndex.find(wtx.hashBlock);
        if (blit != mapBlockIndex.end() && blit->second->IsInMainChain()) {
            // ... which are already in a block
            int nHeight = blit->second->nHeight;
            BOOST_FOREACH(const CTxOut &txout, wtx.vout) {
                // iterate over all their outputs
                ::ExtractAffectedKeys(*this, txout.scriptPubKey, vAffected);
                BOOST_FOREACH(const CKeyID &keyid, vAffected) {
                    // ... and all their affected keys
                    std::map<CKeyID, CBlockIndex*>::iterator rit = mapKeyFirstBlock.find(keyid);
                    if (rit != mapKeyFirstBlock.end() && nHeight < rit->second->nHeight)
                        rit->second = blit->second;
                }
                vAffected.clear();
            }
        }
    }

    // Extract block timestamps for those keys
    for (std::map<CKeyID, CBlockIndex*>::const_iterator it = mapKeyFirstBlock.begin(); it != mapKeyFirstBlock.end(); it++)
        mapKeyBirth[it->first] = it->second->nTime - 7200; // block times can be 2h off
}

void CWallet::SearchNotaryTransactions(uint256 hash, std::vector<std::pair<std::string, int> >& vTxResults)
{
    int blockstogoback = pindexBest->nHeight - 362500;
    std::string matchingCLAMSpeech = "notary " + hash.GetHex();

    const CBlockIndex* pindexFirst = pindexBest;
    for (int i = 0; pindexFirst && i < blockstogoback; i++) {

        CBlock block;
        block.ReadFromDisk(pindexFirst);

        BOOST_FOREACH (const CTransaction& tx, block.vtx)
        {
            if (tx.strCLAMSpeech == matchingCLAMSpeech) {
                vTxResults.push_back( std::make_pair(tx.GetHash().GetHex(), pindexFirst->nHeight) );
            }
        }

        pindexFirst = pindexFirst->pprev;
    }
    return;
}

CClamour* CWallet::GetClamour(std::string pid)
{
    std::map<std::string, CClamour*>::iterator mi = mapClamour.find(pid);
    CClamour *pResult(NULL);
    if (mi == mapClamour.end())
        return pResult;
    pResult = mi->second;
    return pResult;
}

void CWallet::ClearOrphans()
{
    list<uint256> orphans;

    LOCK(cs_wallet);
    for(map<uint256, CWalletTx>::const_iterator it = mapWallet.begin(); it != mapWallet.end(); ++it)
    {
        const CWalletTx *wtx = &(*it).second;
        if((wtx->IsCoinBase() || wtx->IsCoinStake()) && !wtx->IsInMainChain())
        {
            orphans.push_back(wtx->GetHash());
        }
    }

    for(list<uint256>::const_iterator it = orphans.begin(); it != orphans.end(); ++it)
        EraseFromWallet(*it);
}
