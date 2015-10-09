// Copyright (c) 2009-2012 Bitcoin Developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include <iostream>
#include <fstream>

#include "init.h" // for pwalletMain
#include "rpcserver.h"
#include "ui_interface.h"
#include "base58.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/variant/get.hpp>
#include <boost/algorithm/string.hpp>

#include "univalue.h"

using namespace std;

void EnsureWalletIsUnlocked();


namespace bt = boost::posix_time;

// Extended DecodeDumpTime implementation, see this page for details:
// http://stackoverflow.com/questions/3786201/parsing-of-date-time-from-string-boost
const std::locale formats[] = {
    std::locale(std::locale::classic(),new bt::time_input_facet("%Y-%m-%dT%H:%M:%SZ")),
    std::locale(std::locale::classic(),new bt::time_input_facet("%Y-%m-%d %H:%M:%S")),
    std::locale(std::locale::classic(),new bt::time_input_facet("%Y/%m/%d %H:%M:%S")),
    std::locale(std::locale::classic(),new bt::time_input_facet("%d.%m.%Y %H:%M:%S")),
    std::locale(std::locale::classic(),new bt::time_input_facet("%Y-%m-%d"))
};

const size_t formats_n = sizeof(formats)/sizeof(formats[0]);
CWallet* pwalletImport;

std::time_t pt_to_time_t(const bt::ptime& pt)
{
    bt::ptime timet_start(boost::gregorian::date(1970,1,1));
    bt::time_duration diff = pt - timet_start;
    return diff.ticks()/bt::time_duration::rep_type::ticks_per_second;
}

int64_t DecodeDumpTime(const std::string& s)
{
    bt::ptime pt;

    for(size_t i=0; i<formats_n; ++i)
    {
        std::istringstream is(s);
        is.imbue(formats[i]);
        is >> pt;
        if(pt != bt::ptime()) break;
    }

    return pt_to_time_t(pt);
}

std::string static EncodeDumpTime(int64_t nTime) {
    return DateTimeStrFormat("%Y-%m-%dT%H:%M:%SZ", nTime);
}

std::string static EncodeDumpString(const std::string &str) {
    std::stringstream ret;
    BOOST_FOREACH(unsigned char c, str) {
        if (c <= 32 || c >= 128 || c == '%') {
            ret << '%' << HexStr(&c, &c + 1);
        } else {
            ret << c;
        }
    }
    return ret.str();
}

std::string DecodeDumpString(const std::string &str) {
    std::stringstream ret;
    for (unsigned int pos = 0; pos < str.length(); pos++) {
        unsigned char c = str[pos];
        if (c == '%' && pos+2 < str.length()) {
            c = (((str[pos+1]>>6)*9+((str[pos+1]-'0')&15)) << 4) | 
                ((str[pos+2]>>6)*9+((str[pos+2]-'0')&15));
            pos += 2;
        }
        ret << c;
    }
    return ret.str();
}

class CTxDump
{
public:
    CBlockIndex *pindex;
    int64_t nValue;
    bool fSpent;
    CWalletTx* ptx;
    int nOut;
    CTxDump(CWalletTx* ptx = NULL, int nOut = -1)
    {
        pindex = NULL;
        nValue = 0;
        fSpent = false;
        this->ptx = ptx;
        this->nOut = nOut;
    }
};

UniValue importprivkey(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 3)
        throw runtime_error(
            "importprivkey <clamprivkey> [label] [rescan=true]\n"
            "Adds a private key (as returned by dumpprivkey) to your wallet.\n"
	    "Accepts CLAM, BTC, LTC and DOGE private keys\n");

    string strSecret = params[0].get_str();
    string strLabel = "";
    if (params.size() > 1)
        strLabel = params[1].get_str();

    // Whether to perform rescan after import
    bool fRescan = true;
    if (params.size() > 2)
        fRescan = params[2].get_bool();

    CBitcoinSecret vchSecret;
    bool fGood = vchSecret.SetString(strSecret);

    if (!fGood) throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid private key");
    if (fWalletUnlockStakingOnly)
        throw JSONRPCError(RPC_WALLET_UNLOCK_NEEDED, "Wallet is unlocked for staking only.");

    CKey key = vchSecret.GetKey();
    CPubKey pubkey = key.GetPubKey();
    CKeyID vchAddress = pubkey.GetID();
    {
        LOCK2(cs_main, pwalletMain->cs_wallet);

        pwalletMain->MarkDirty();
        pwalletMain->SetAddressBookName(vchAddress, strLabel);

        // Don't throw error in case a key is already there
        if (pwalletMain->HaveKey(vchAddress))
            return NullUniValue;

        pwalletMain->mapKeyMetadata[vchAddress].nCreateTime = 1;

        if (!pwalletMain->AddKeyPubKey(key, pubkey))
            throw JSONRPCError(RPC_WALLET_ERROR, "Error adding key to wallet");

        // whenever a key is imported, we need to scan the whole chain
        pwalletMain->nTimeFirstKey = 1; // 0 would be considered 'no value'

        if (fRescan) {
            pwalletMain->ScanForWalletTransactions(pindexGenesisBlock, true);
            pwalletMain->ReacceptWalletTransactions();
        }
    }

    return NullUniValue;
}

UniValue deleteprivkey(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "deleteprivkey <clamaddress>\n"
            "Irreversibly deletes the private key corresponding to the given address from your wallet.\n"
	    "Accepts CLAM addresses\n");

    EnsureWalletIsUnlocked();

    string strAddress = params[0].get_str();
    CKeyID vchAddress;
    CPubKey vchPubKey;

    if (!CBitcoinAddress(strAddress).GetKeyID(vchAddress))
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid address");

    if (fWalletUnlockStakingOnly)
        throw JSONRPCError(RPC_WALLET_UNLOCK_NEEDED, "Wallet is unlocked for staking only.");

    {
        LOCK2(cs_main, pwalletMain->cs_wallet);

        if (!pwalletMain->HaveKey(vchAddress))
            throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "We don't have the private key for that address");

        if (!pwalletMain->GetPubKey(vchAddress, vchPubKey))
            throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Can't find pubkey for that address");
            
        pwalletMain->MarkDirty();
        pwalletMain->DelAddressBookName(vchAddress);
        pwalletMain->mapKeyMetadata.erase(vchAddress);
        pwalletMain->RemovePubKey(vchPubKey);
        pwalletMain->TidyWalletTransactions();
    }

    return NullUniValue;
}

UniValue importwalletdump(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "importwalletdump <filename>\n"
            "Imports keys from a wallet dump file (see dumpwallet).");

    EnsureWalletIsUnlocked();

    ifstream file;
    file.open(params[0].get_str().c_str());
    if (!file.is_open())
        throw JSONRPCError(RPC_INVALID_PARAMETER, "Cannot open wallet dump file");

    int64_t nTimeBegin = pindexBest->nTime;

    bool fGood = true;

    while (file.good()) {
        std::string line;
        std::getline(file, line);
        if (line.empty() || line[0] == '#')
            continue;

        std::vector<std::string> vstr;
        boost::split(vstr, line, boost::is_any_of(" "));
        if (vstr.size() < 2)
            continue;
        CBitcoinSecret vchSecret;
        if (!vchSecret.SetString(vstr[0]))
            continue;
        CKey key = vchSecret.GetKey();
        CPubKey pubkey = key.GetPubKey();
        CKeyID keyid = pubkey.GetID();
        if (pwalletMain->HaveKey(keyid)) {
            LogPrintf("Skipping import of %s (key already present)\n", CBitcoinAddress(keyid).ToString());
            continue;
        }
        int64_t nTime = DecodeDumpTime(vstr[1]);
        std::string strLabel;
        bool fLabel = true;
        for (unsigned int nStr = 2; nStr < vstr.size(); nStr++) {
            if (boost::algorithm::starts_with(vstr[nStr], "#"))
                break;
            if (vstr[nStr] == "change=1")
                fLabel = false;
            if (vstr[nStr] == "reserve=1")
                fLabel = false;
            if (boost::algorithm::starts_with(vstr[nStr], "label=")) {
                strLabel = DecodeDumpString(vstr[nStr].substr(6));
                fLabel = true;
            }
        }
        LogPrintf("Importing %s...\n", CBitcoinAddress(keyid).ToString());
        if (!pwalletMain->AddKey(key)) {
            fGood = false;
            continue;
        }
        pwalletMain->mapKeyMetadata[keyid].nCreateTime = nTime;
        if (fLabel)
            pwalletMain->SetAddressBookName(keyid, strLabel);
        nTimeBegin = std::min(nTimeBegin, nTime);
    }
    file.close();

    CBlockIndex *pindex = pindexBest;
    while (pindex && pindex->pprev && pindex->nTime > nTimeBegin - 7200)
        pindex = pindex->pprev;

    if (!pwalletMain->nTimeFirstKey || nTimeBegin < pwalletMain->nTimeFirstKey)
        pwalletMain->nTimeFirstKey = nTimeBegin;

    LogPrintf("Rescanning last %i blocks\n", pindexBest->nHeight - pindex->nHeight + 1);
    pwalletMain->ScanForWalletTransactions(pindex);
    pwalletMain->ReacceptWalletTransactions();
    pwalletMain->MarkDirty();

    if (!fGood)
        throw JSONRPCError(RPC_WALLET_ERROR, "Error adding some keys to wallet");

    return NullUniValue;
}


UniValue dumpprivkey(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "dumpprivkey <clamaddress>\n"
            "Reveals the private key corresponding to <clamaddress>.");

    EnsureWalletIsUnlocked();

    string strAddress = params[0].get_str();
    CBitcoinAddress address;
    if (!address.SetString(strAddress))
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid Clam address");
    if (fWalletUnlockStakingOnly)
        throw JSONRPCError(RPC_WALLET_UNLOCK_NEEDED, "Wallet is unlocked for staking only.");
    CKeyID keyID;
    if (!address.GetKeyID(keyID))
        throw JSONRPCError(RPC_TYPE_ERROR, "Address does not refer to a key");
    CKey vchSecret;
    if (!pwalletMain->GetKey(keyID, vchSecret))
        throw JSONRPCError(RPC_WALLET_ERROR, "Private key for address " + strAddress + " is not known");
    return CBitcoinSecret(vchSecret).ToString();
}

UniValue dumpwallet(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "dumpwallet <filename>\n"
            "Dumps all wallet keys in a human-readable format.");

    EnsureWalletIsUnlocked();

    ofstream file;
    file.open(params[0].get_str().c_str());
    if (!file.is_open())
        throw JSONRPCError(RPC_INVALID_PARAMETER, "Cannot open wallet dump file");

    std::map<CKeyID, int64_t> mapKeyBirth;

    std::set<CKeyID> setKeyPool;

    pwalletMain->GetKeyBirthTimes(mapKeyBirth);

    pwalletMain->GetAllReserveKeys(setKeyPool);

    // sort time/key pairs
    std::vector<std::pair<int64_t, CKeyID> > vKeyBirth;
    for (std::map<CKeyID, int64_t>::const_iterator it = mapKeyBirth.begin(); it != mapKeyBirth.end(); it++) {
        vKeyBirth.push_back(std::make_pair(it->second, it->first));
    }
    mapKeyBirth.clear();
    std::sort(vKeyBirth.begin(), vKeyBirth.end());

    // produce output
    file << strprintf("# Wallet dump created by Clam %s (%s)\n", CLIENT_BUILD, CLIENT_DATE);
    file << strprintf("# * Created on %s\n", EncodeDumpTime(GetTime()));
    file << strprintf("# * Best block at time of backup was %i (%s),\n", nBestHeight, hashBestChain.ToString());
    file << strprintf("#   mined on %s\n", EncodeDumpTime(pindexBest->nTime));
    file << "\n";
    for (std::vector<std::pair<int64_t, CKeyID> >::const_iterator it = vKeyBirth.begin(); it != vKeyBirth.end(); it++) {
        const CKeyID &keyid = it->second;
        std::string strTime = EncodeDumpTime(it->first);
        std::string strAddr = CBitcoinAddress(keyid).ToString();

        CKey key;
        if (pwalletMain->GetKey(keyid, key)) {
            if (pwalletMain->mapAddressBook.count(keyid)) {
                file << strprintf("%s %s label=%s # addr=%s\n", CBitcoinSecret(key).ToString(), strTime, EncodeDumpString(pwalletMain->mapAddressBook[keyid]), strAddr);
            } else if (setKeyPool.count(keyid)) {
                file << strprintf("%s %s reserve=1 # addr=%s\n", CBitcoinSecret(key).ToString(), strTime, strAddr);
            } else {
                file << strprintf("%s %s change=1 # addr=%s\n", CBitcoinSecret(key).ToString(), strTime, strAddr);
            }
        }
    }
    file << "\n";
    file << "# End of dump\n";
    file.close();
    return NullUniValue;
}

UniValue importwallet(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 3)
        throw runtime_error(
            "importwallet <file> [walletpassword] [rescan=true]\n"
            "Import wallet.dat from BTC/LTC/DOGE/CLAM \n"
            "Password is only required if wallet is encrypted\n"
        );

    // Whether to perform rescan after import
    bool fRescan = true;
    if (params.size() > 2)
        fRescan = params[2].get_bool();

    CBlockIndex *pindexRescan = pindexGenesisBlock;
    EnsureWalletIsUnlocked();

    bool fFirstRun = false;

    pwalletImport = new CWallet(params[0].get_str().c_str());
    DBErrors nLoadWalletRet = pwalletImport->LoadWalletImport(fFirstRun);

    std::ostringstream strErrors;
    if (nLoadWalletRet != DB_LOAD_OK)
        {
            if (nLoadWalletRet == DB_CORRUPT)
                throw JSONRPCError(RPC_WALLET_ERROR, "Error loading wallet.dat: Wallet corrupted");
            else if (nLoadWalletRet == DB_LOAD_FAIL)
            {
                throw JSONRPCError(RPC_WALLET_ERROR, "Wallet failed to load");
            }
            else
                LogPrintf("Non-fatal errors loading wallet file\n");
    }

    // Handle encrypted wallets. Wallets first need to be unlocked before the keys
    // can be added into your clam wallet. 
    if (pwalletImport->IsCrypted() && pwalletImport->IsLocked()) {
        bool fGotWalletPass = true;
        if (params.size() < 2)
            fGotWalletPass = false;
        else
        {
            // TODO: get rid of this .c_str() by implementing SecureString::operator=(std::string)
            // Alternately, find a way to make params[0] mlock()'d to begin with.
            SecureString strWalletPass;
            strWalletPass.reserve(100);
            strWalletPass = params[1].get_str().c_str();
            if (strWalletPass.length() > 0)
            {
                if (!pwalletImport->Unlock(strWalletPass))
                    throw JSONRPCError(RPC_WALLET_PASSPHRASE_INCORRECT, "Error: The wallet passphrase entered was incorrect for the wallet you are attempting to import.");
            } else
                fGotWalletPass = false;
        }

        if (!fGotWalletPass)
            throw runtime_error(
                "importwallet <file> <walletpassword>\n"
                "Import encrypted wallet from BTC/LTC/DOGE \n\n"
                "You are attempting to import an encrypted wallet\n"
                "The passphrase must be entered to import the wallet\n"
                );
    }

    int nImported = 0;
    int nSkipped = 0;

    {
        LOCK2(cs_main, pwalletMain->cs_wallet);
        LOCK(pwalletImport->cs_wallet);

        std::set<CKeyID> setKeys;
        pwalletImport->GetKeys(setKeys);

        BOOST_FOREACH(const CKeyID &keyid, setKeys) {
            
            int64_t nTime = GetTime();

            std::string strAddr = CBitcoinAddress(keyid).ToString();
            std::string strLabel = "importwallet";

            CKey key;
            if (pwalletImport->GetKey(keyid, key)) {

                if (pwalletMain->HaveKey(keyid)) {
                    if (fDebug) 
                        LogPrintf("Skipping import of %s (key already present)\n", strAddr);
                    nSkipped++;
                    continue;
                }
		
                if (fDebug) 
                    LogPrintf("Importing %s...\n", strAddr);

                pwalletMain->AddKey(key);
                pwalletMain->SetAddressBookName(keyid, strLabel);
                pwalletMain->mapKeyMetadata[keyid].nCreateTime = nTime;
                nImported++;
            }
        }
    }
    
    // Clean up unregistered wallet
    UnregisterWallet(pwalletImport);
    delete pwalletImport;

    LogPrintf("walletimport imported %d and skipped %d key(s)\n", nImported, nSkipped);
    pwalletMain->MarkDirty();

    if (nImported)
        if (fRescan)
        {
            LogPrintf("Searching last %i blocks (from block %i) for dug Clams...\n", pindexBest->nHeight - pindexRescan->nHeight, pindexRescan->nHeight);
            pwalletMain->ScanForWalletTransactions(pindexRescan, true);
            pwalletMain->ReacceptWalletTransactions();
            LogPrintf("Rescan complete\n");
        }
        else
            LogPrintf("Not rescanning because user requested that it should be skipped\n");
    else
        LogPrintf("Not rescanning because no new keys were imported\n");

    return NullUniValue;
}
