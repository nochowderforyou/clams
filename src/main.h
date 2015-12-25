// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.
#ifndef BITCOIN_MAIN_H
#define BITCOIN_MAIN_H

#if defined(HAVE_CONFIG_H)
#include "bitcoin-config.h"
#endif

#include "core.h"
#include "bignum.h"
#include "sync.h"
#include "net.h"
#include "script.h"
#include "scrypt.h"

#include <limits>
#include <list>

class CBlock;
class CBlockIndex;
class CInv;
class CKeyItem;
class CNode;
class CReserveKey;
class CWallet;

struct CBlockIndexWorkComparator;

static const int LAST_POW_BLOCK = 10000;
static const int LAST_TESTNET_POW_BLOCK = 300;

/** The maximum allowed size for a serialized block, in bytes (network rule) */
static const unsigned int MAX_BLOCK_SIZE = 1000000;
/** The maximum size for mined blocks */
static const unsigned int MAX_BLOCK_SIZE_GEN = MAX_BLOCK_SIZE/2;
/** The maximum size for transactions we're willing to relay/mine **/
static const unsigned int MAX_STANDARD_TX_SIZE = MAX_BLOCK_SIZE_GEN/5;
/** The maximum allowed number of signature check operations in a block (network rule) */
static const unsigned int MAX_BLOCK_SIGOPS = MAX_BLOCK_SIZE/50;
/** The maximum number of orphan transactions kept in memory */
static const unsigned int MAX_ORPHAN_TRANSACTIONS = MAX_BLOCK_SIZE/100;
/** The maximum number of entries in an 'inv' protocol message */
static const unsigned int MAX_INV_SZ = 50000;
/** The maximum size of a blk?????.dat file (since 0.8) */
static const unsigned int MAX_BLOCKFILE_SIZE = 0x8000000; // 128 MiB
/** The pre-allocation chunk size for blk?????.dat files (since 0.8) */
static const unsigned int BLOCKFILE_CHUNK_SIZE = 0x1000000; // 16 MiB
/** The pre-allocation chunk size for rev?????.dat files (since 0.8) */
static const unsigned int UNDOFILE_CHUNK_SIZE = 0x100000; // 1 MiB
/** The maximum number of sigops we're willing to relay/mine in a single tx */
static const unsigned int MAX_TX_SIGOPS = MAX_BLOCK_SIGOPS/5;
/** Default for -maxorphanblocks, maximum number of orphan blocks kept in memory */
static const unsigned int DEFAULT_MAX_ORPHAN_BLOCKS = 750;
/** Fees smaller than this (in satoshi) are considered zero fee (for transaction creation) */
static const int64_t MIN_TX_FEE = 10000;
/** Fees smaller than this (in satoshi) are considered zero fee (for relaying) */
static const int64_t MIN_RELAY_TX_FEE = MIN_TX_FEE;
static const int64_t MAX_MONEY = std::numeric_limits<int64_t>::max();
static const int64_t COIN_YEAR_REWARD = 1 * CENT; // 1% per year

inline bool MoneyRange(int64_t nValue) { return (nValue >= 0 && nValue <= MAX_MONEY); }
/** Threshold for nLockTime: below this value it is interpreted as block number, otherwise as UNIX timestamp. */
static const unsigned int LOCKTIME_THRESHOLD = 500000000; // Tue Nov  5 00:53:20 1985 UTC
/** Maximum number of script-checking threads allowed */
static const int MAX_SCRIPTCHECK_THREADS = 16;
/** Fake height value used in CCoins to signify they are only in the memory pool (since 0.8) */
static const unsigned int MEMPOOL_HEIGHT = 0x7FFFFFFF;

/** Number bytes a transaction comment is allowed to have */
static const unsigned int MAX_TX_COMMENT_LEN = 140; // 128 bytes + little extra

static const uint256 hashGenesisBlock("0x00000c3ce6b3d823a35224a39798eca9ad889966aeb5a9da7b960ffb9869db35");
static const uint256 hashGenesisBlockTestNet("0x00001924120e93f445dd4adb9d90e0020350b8c6c2b08e1a4950372a37f8bcc8");

static const uint256 hashHighBlock ("0xdb61f591d7fb40afa08476d6492e81a06edddf332d7027968ac130db95c07cb7");
static const int HIGH_BLOCK_INDEX = 275000;


inline bool IsProtocolV2(int nHeight) { return TestNet() || nHeight > 203500; }

inline int64_t PastDrift(int64_t nTime, int nHeight)   { return IsProtocolV2(nHeight) ? nTime      : nTime - 10 * 60; }

inline int64_t FutureDriftV1(int64_t nTime) { return nTime + 10 * 60; }
inline int64_t FutureDriftV2(int64_t nTime) { return nTime + 15; }
inline int64_t FutureDrift(int64_t nTime, int nHeight) { return IsProtocolV2(nHeight) ? FutureDriftV2(nTime) : FutureDriftV1(nTime); }

static const int64 MAX_TIME_SINCE_BEST_BLOCK = 10; // how many seconds to wait before sending next PushGetBlocks()
extern CScript COINBASE_FLAGS;
extern CCriticalSection cs_main;
extern std::map<uint256, CBlockIndex*> mapBlockIndex;
extern std::set<CBlockIndex*, CBlockIndexWorkComparator> setBlockIndexValid;
extern std::set<std::pair<COutPoint, unsigned int> > setStakeSeen;
extern CBlockIndex* pindexGenesisBlock;
extern unsigned int nTargetSpacing;
extern unsigned int nTargetStakeSpacing;
extern unsigned int nStakeMinAge;
extern unsigned int nStakeMaxAge;
extern unsigned int nNodeLifespan;
extern int nCoinbaseMaturity;
extern int nBestHeight;
extern uint256 nBestChainWork;
extern uint256 nBestInvalidWork;
extern uint256 hashBestChain;
extern CBlockIndex* pindexBest;
extern unsigned int nTransactionsUpdated;
extern uint64_t nLastBlockTx;
extern uint64_t nLastBlockSize;
extern int64_t nLastCoinStakeSearchInterval;
extern const std::string strMessageMagic;
extern int64_t nTimeBestReceived;
extern bool fImporting;
extern bool fReindex;
extern bool fBenchmark;
extern int nScriptCheckThreads;
struct COrphanBlock;
extern std::map<uint256, COrphanBlock*> mapOrphanBlocks;
extern std::map<uint256, CBlock*> mapDuplicateStakeBlocks;
extern bool fTxIndex;
extern unsigned int nCoinCacheSize;
extern bool fHaveGUI;

extern int nBlocksToIgnore;

// Settings
extern bool fUseFastIndex;
extern unsigned int nDerivationMethodIndex;

extern bool fMinimizeCoinAge;
extern bool fCreditStakesToAccounts;

// Minimum disk space required - used in CheckDiskSpace()
static const uint64_t nMinDiskSpace = 52428800;

class CReserveKey;
class CCoinsDB;
class CBlockTreeDB;
struct CDiskBlockPos;
class CCoins;
class CTxUndo;
class CCoinsView;
class CCoinsViewCache;
class CScriptCheck;
class CValidationState;
class CWalletInterface;

struct CBlockTemplate;

/** Register a wallet to receive updates from core */
void RegisterWallet(CWalletInterface* pwalletIn);
/** Unregister a wallet from core */
void UnregisterWallet(CWalletInterface* pwalletIn);
/** Unregister all wallets from core */
void UnregisterAllWallets();
/** Push an updated transaction to all registered wallets */
void SyncWithWallets(const CTransaction& tx, const CBlock* pblock = NULL, bool fConnect = true);
/** Inform all registered wallets of a change in stake reward */
void StakeToWallets(const CScript& script, int64_t nStakeReward, bool fConnect = true);
/** Ask wallets to resend their transactions */
void ResendWalletTransactions(bool fForce = false);

/** Register with a network node to receive its signals */
void RegisterNodeSignals(CNodeSignals& nodeSignals);
/** Unregister a network node */
void UnregisterNodeSignals(CNodeSignals& nodeSignals);

void PushGetBlocks(CNode* pnode, CBlockIndex* pindexBegin, uint256 hashEnd);

bool ProcessBlock(CValidationState &state, CNode* pfrom, CBlock* pblock, CDiskBlockPos *dbp = NULL);
bool CheckDiskSpace(uint64_t nAdditionalBytes=0);

/** Open a block file (blk?????.dat) */
FILE* OpenBlockFile(const CDiskBlockPos &pos, bool fReadOnly = false);
/** Open an undo file (rev?????.dat) */
FILE* OpenUndoFile(const CDiskBlockPos &pos, bool fReadOnly = false);
/** Import blocks from an external file */
bool LoadExternalBlockFile(FILE* fileIn, CDiskBlockPos *dbp = NULL);
/** Initialize a new block tree database + block data on disk */
bool InitBlockIndex();
/** Load the block tree and coins database from disk */
bool LoadBlockIndex();
void UnloadBlockIndex();
/** Verify consistency of the block and coin databases */
bool VerifyDB(int nCheckLevel, int nCheckDepth);
/** Print the loaded block tree */
void PrintBlockTree();
CBlockIndex* FindBlockByHeight(int nHeight);
bool ProcessMessages(CNode* pfrom);
bool SendMessages(CNode* pto, bool fSendTrickle);

bool CheckProofOfWork(uint256 hash, unsigned int nBits);
unsigned int GetNextTargetRequired(const CBlockIndex* pindexLast, bool fProofOfStake);
int64_t GetProofOfWorkReward(int nHeight, int64_t nFees);
int64_t GetProofOfStakeReward(const CBlockIndex* pindex, int64_t nCoinAge, int64_t nFees);
unsigned int ComputeMinWork(unsigned int nBase, int64_t nTime);
unsigned int ComputeMinStake(unsigned int nBase, int64_t nTime, unsigned int nBlockTime);
bool IsInitialBlockDownload();
std::string GetWarnings(std::string strFor);
/** Retrieve a transaction (from memory pool, or from disk, if possible) */
bool GetTransaction(const uint256 &hash, CTransaction &tx, uint256 &hashBlock, bool fAllowSlow = false);
/** Connect/disconnect blocks until pindexNew is the new tip of the active block chain */
bool SetBestChain(CValidationState &state, CBlockIndex* pindexNew);
/** Find the best known block, and make it the tip of the block chain */
bool ConnectBestBlock(CValidationState &state);
/** Create a new block index entry for a given block hash */
CBlockIndex * InsertBlockIndex(uint256 hash);
uint256 WantedByOrphan(const COrphanBlock* pblockOrphan);
const CBlockIndex* GetLastBlockIndex(const CBlockIndex* pindex, bool fProofOfStake);
void ThreadStakeMiner(CWallet *pwallet);
/** Abort with a message */
bool AbortNode(const std::string &msg);
/** Verify a signature */
bool VerifySignature(const CTransaction& txFrom, const CTransaction& txTo, unsigned int nIn, int nHashType, int nBlockHeight, int nBlockTime);






struct CDiskBlockPos
{
    int nFile;
    unsigned int nPos;

    IMPLEMENT_SERIALIZE(
        READWRITE(VARINT(nFile));
        READWRITE(VARINT(nPos));
    )

    CDiskBlockPos() {
        SetNull();
    }

    CDiskBlockPos(int nFileIn, unsigned int nPosIn) {
        nFile = nFileIn;
        nPos = nPosIn;
    }

    friend bool operator==(const CDiskBlockPos &a, const CDiskBlockPos &b) {
        return (a.nFile == b.nFile && a.nPos == b.nPos);
    }

    friend bool operator!=(const CDiskBlockPos &a, const CDiskBlockPos &b) {
        return !(a == b);
    }

    void SetNull() { nFile = -1; nPos = 0; }
    bool IsNull() const { return (nFile == -1); }
};

struct CDiskTxPos : public CDiskBlockPos
{
    unsigned int nTxOffset; // after header
     unsigned int nTxPos;

    IMPLEMENT_SERIALIZE(
        READWRITE(*(CDiskBlockPos*)this);
        READWRITE(VARINT(nTxOffset));
        READWRITE(VARINT(nTxPos));
    )

    CDiskTxPos(const CDiskBlockPos &blockIn, unsigned int nTxOffsetIn, unsigned int nTxPos) : CDiskBlockPos(blockIn.nFile, blockIn.nPos), nTxOffset(nTxOffsetIn), nTxPos(nTxPos) {
    }

    CDiskTxPos() {
        SetNull();
    }

    void SetNull() {
        CDiskBlockPos::SetNull();
        nTxOffset = 0;
        nTxPos =0;
    }
};


bool GetTransaction(const uint256 &hash, CTransaction &tx, uint256 &hashBlock, CDiskTxPos &postx, bool fAllowSlow = false);



enum GetMinFee_mode
{
    GMF_BLOCK,
    GMF_RELAY,
    GMF_SEND,
};

int64_t GetMinFee(const CTransaction& tx, unsigned int nBlockSize = 1, enum GetMinFee_mode mode = GMF_BLOCK, unsigned int nBytes = 0);

/** The basic transaction that is broadcasted on the network and contained in
 * blocks.  A transaction can contain multiple inputs and outputs.
 */
class CTransaction
{
public:
    static int64 nMinTxFee;
    static int64 nMinRelayTxFee;
    static const int CURRENT_VERSION=2;
    static const int LEGACY_VERSION_1=1;
    int nVersion;
    unsigned int nTime;
    std::vector<CTxIn> vin;
    std::vector<CTxOut> vout;
    unsigned int nLockTime;
    std::string strCLAMSpeech;
    mutable uint256 hash;

    // Denial-of-service detection:
    mutable int nDoS;
    bool DoS(int nDoSIn, bool fIn) const { nDoS += nDoSIn; return fIn; }

    CTransaction()
    {
        SetNull();
    }

    IMPLEMENT_SERIALIZE
    (
        READWRITE(this->nVersion);
        nVersion = this->nVersion;
        READWRITE(nTime);
        READWRITE(vin);
        READWRITE(vout);
        READWRITE(nLockTime);
        if (this->nVersion > LEGACY_VERSION_1)
        {
          READWRITE(strCLAMSpeech);
        }
    )

    void SetNull()
    {
        //nVersion = CTransaction::CURRENT_VERSION;
        nVersion = IsProtocolV2(nBestHeight) ? CTransaction::CURRENT_VERSION : LEGACY_VERSION_1;
        nTime = GetAdjustedTime();
        vin.clear();
        vout.clear();
        nLockTime = 0;
        nDoS = 0;  // Denial-of-service prevention
        strCLAMSpeech.clear();
        hash = 0;
    }

    bool IsNull() const
    {
        return (vin.empty() && vout.empty());
    }

    uint256 GetHash() const
    {
        return SerializeHash(*this);
    }

    bool IsNewerThan(const CTransaction& old) const
    {
        if (vin.size() != old.vin.size())
            return false;
        for (unsigned int i = 0; i < vin.size(); i++)
            if (vin[i].prevout != old.vin[i].prevout)
                return false;

        bool fNewer = false;
        unsigned int nLowest = std::numeric_limits<unsigned int>::max();
        for (unsigned int i = 0; i < vin.size(); i++)
        {
            if (vin[i].nSequence != old.vin[i].nSequence)
            {
                if (vin[i].nSequence <= nLowest)
                {
                    fNewer = false;
                    nLowest = vin[i].nSequence;
                }
                if (old.vin[i].nSequence < nLowest)
                {
                    fNewer = true;
                    nLowest = old.vin[i].nSequence;
                }
            }
        }
        return fNewer;
    }

    bool IsCoinBase() const
    {
        return (vin.size() == 1 && vin[0].prevout.IsNull() && vout.size() >= 1);
    }

    bool IsCoinStake() const
    {
        // ppcoin: the coin stake transaction is marked with the first output empty
        return (vin.size() > 0 && (!vin[0].prevout.IsNull()) && vout.size() >= 2 && vout[0].IsEmpty());
    }

    bool IsCreateClamour(std::string& strHash, std::string& strURL) const;
    /** Check for standard transaction types
    @return True if all outputs (scriptPubKeys) use only standard transaction forms
    */
    bool IsStandard(std::string& strReason) const;
    bool IsStandard() const
    {
        std::string strReason;
        return IsStandard(strReason);
    }

    /** Check for standard transaction types
    @param[in] mapInputs    Map of previous transactions that have outputs we're spending
    @return True if all inputs (scriptSigs) use only standard transaction forms
    */
    bool AreInputsStandard(CCoinsViewCache& mapInputs) const;

    /** Count ECDSA signature operations the old-fashioned (pre-0.6) way
        @return number of sigops this transaction's outputs will produce when spent
    */
    unsigned int GetLegacySigOpCount() const;

    /** Count ECDSA signature operations in pay-to-script-hash inputs.

        @param[in] mapInputs    Map of previous transactions that have outputs we're spending
        @return maximum number of sigops required to validate this transaction's inputs
     */
    unsigned int GetP2SHSigOpCount(CCoinsViewCache& mapInputs) const;

    /** Amount of bitcoins spent by this transaction.
        @return sum of all outputs (note: does not include fees)
     */
    int64_t GetValueOut() const
    {
        int64_t nValueOut = 0;
        BOOST_FOREACH(const CTxOut& txout, vout)
        {
            nValueOut += txout.nValue;
            if (!MoneyRange(txout.nValue) || !MoneyRange(nValueOut))
                throw std::runtime_error("CTransaction::GetValueOut() : value out of range");
        }
        return nValueOut;
    }

    /** Amount of bitcoins coming in to this transaction
        Note that lightweight clients may not know anything besides the hash of previous transactions,
        so may not be able to calculate this.

        @param[in] mapInputs	Map of previous transactions that have outputs we're spending
        @return	Sum of value of all inputs (scriptSigs)
        @see CTransaction::FetchInputs
     */
    int64_t GetValueIn(CCoinsViewCache& mapInputs) const;

    // Apply the effects of this transaction on the UTXO set represented by view
    void UpdateCoins(const CTransaction& tx, CValidationState &state, CCoinsViewCache &inputs, CTxUndo &txundo, int nHeight, const uint256 &txhash);

    int64_t GetMinFee(unsigned int nBlockSize=1, bool fAllowFree=true, enum GetMinFee_mode mode=GMF_BLOCK) const;

    friend bool operator==(const CTransaction& a, const CTransaction& b)
    {
        return (a.nVersion  == b.nVersion &&
                a.nTime     == b.nTime &&
                a.vin       == b.vin &&
                a.vout      == b.vout &&
                a.nLockTime == b.nLockTime);
    }

    friend bool operator!=(const CTransaction& a, const CTransaction& b)
    {
        return !(a == b);
    }

    std::string ToString() const
    {
        std::string str;
        str += IsCoinBase()? "Coinbase" : (IsCoinStake()? "Coinstake" : "CTransaction");
        str += strprintf("(hash=%s, nTime=%d, ver=%d, vin.size=%"PRIszu", vout.size=%"PRIszu", nLockTime=%d, strCLAMSpeech=%s)\n",
            GetHash().ToString().substr(0,10).c_str(),
            nTime,
            nVersion,
            vin.size(),
            vout.size(),
            nLockTime,
            strCLAMSpeech.substr(0,30).c_str());
        for (unsigned int i = 0; i < vin.size(); i++)
            str += "    " + vin[i].ToString() + "\n";
        for (unsigned int i = 0; i < vout.size(); i++)
            str += "    " + vout[i].ToString() + "\n";
        return str;
    }

    void print() const
    {
        LogPrintf("%s", ToString().c_str());
    }

    // Check whether all prevouts of this transaction are present in the UTXO set represented by view
    bool HaveInputs(CCoinsViewCache &view, int64_t* nDigs = 0) const;

    // Check whether all inputs of this transaction are valid (no double spends, scripts & sigs, amounts)
    // This does not modify the UTXO set. If pvChecks is not NULL, script checks are pushed onto it
    // instead of being performed inline.
    bool CheckInputs(CValidationState &state, CCoinsViewCache &view, bool fScriptChecks = true,
                     unsigned int flags = SCRIPT_VERIFY_NOCACHE | SCRIPT_VERIFY_LOCKTIME) const;

    // Apply the effects of this transaction on the UTXO set represented by view
    void UpdateCoins(CValidationState &state, CCoinsViewCache &view, CTxUndo &txundo, int nHeight, const uint256 &txhash) const;

    // Context-independent validity checks
    bool CheckTransaction(CValidationState &state) const;

    // Try to accept this transaction into the memory pool
    bool AcceptToMemoryPool(CValidationState &state, bool fCheckInputs=true, bool fLimitFree = true, bool* pfMissingInputs=NULL, bool fRejectInsaneFee = false);

    // ppcoin: get transaction coin age
    bool GetCoinAge(uint64_t& nCoinAge) const;

protected:
    static const CTxOut &GetOutputFor(const CTxIn& input, CCoinsViewCache& mapInputs);
};


bool IsFinalTx(const CTransaction &tx, int nBlockHeight = 0, int64_t nBlockTime = 0);
bool IsFinalLockTime(unsigned int nLockTime, int nBlockHeight = 0, int64_t nBlockTime = 0);



/** wrapper for CTxOut that provides a more compact serialization */
class CTxOutCompressor
{
private:
    CTxOut &txout;
public:
    CTxOutCompressor(CTxOut &txoutIn) : txout(txoutIn) { }

    IMPLEMENT_SERIALIZE(
        READWRITE(VARINT(txout.nValue));
        CScriptCompressor cscript(REF(txout.scriptPubKey));
        READWRITE(cscript);
    )
};

/** Undo information for a CTxIn
 *
 *  Contains the prevout's CTxOut being spent, and if this was the
 *  last output of the affected transaction, its metadata as well
 *  (coinbase or not, height, transaction version)
 */
class CTxInUndo
{
public:
    CTxOut txout;         // the txout data before being spent
    bool fCoinBase;       // if the outpoint was the last unspent: whether it belonged to a coinbase
    unsigned int nHeight; // if the outpoint was the last unspent: its height
    int nVersion;         // if the outpoint was the last unspent: its version

    CTxInUndo() : txout(), fCoinBase(false), nHeight(0), nVersion(0) {}
    CTxInUndo(const CTxOut &txoutIn, bool fCoinBaseIn = false, unsigned int nHeightIn = 0, int nVersionIn = 0) : txout(txoutIn), fCoinBase(fCoinBaseIn), nHeight(nHeightIn), nVersion(nVersionIn) { }

    unsigned int GetSerializeSize(int nType, int nVersion) const {
        return ::GetSerializeSize(VARINT(nHeight*2+(fCoinBase ? 1 : 0)), nType, nVersion) +
               (nHeight > 0 ? ::GetSerializeSize(VARINT(this->nVersion), nType, nVersion) : 0) +
               ::GetSerializeSize(CTxOutCompressor(REF(txout)), nType, nVersion);
    }

    template<typename Stream>
    void Serialize(Stream &s, int nType, int nVersion) const {
        ::Serialize(s, VARINT(nHeight*2+(fCoinBase ? 1 : 0)), nType, nVersion);
        if (nHeight > 0)
            ::Serialize(s, VARINT(this->nVersion), nType, nVersion);
        ::Serialize(s, CTxOutCompressor(REF(txout)), nType, nVersion);
    }

    template<typename Stream>
    void Unserialize(Stream &s, int nType, int nVersion) {
        unsigned int nCode = 0;
        ::Unserialize(s, VARINT(nCode), nType, nVersion);
        nHeight = nCode / 2;
        fCoinBase = nCode & 1;
        if (nHeight > 0)
            ::Unserialize(s, VARINT(this->nVersion), nType, nVersion);
        ::Unserialize(s, REF(CTxOutCompressor(REF(txout))), nType, nVersion);
    }
};

/** Undo information for a CTransaction */
class CTxUndo
{
public:
    // undo information for all txins
    std::vector<CTxInUndo> vprevout;

    IMPLEMENT_SERIALIZE(
        READWRITE(vprevout);
    )
};

/** Undo information for a CBlock */
class CBlockUndo
{
public:
    std::vector<CTxUndo> vtxundo; // for all but the coinbase

    IMPLEMENT_SERIALIZE(
        READWRITE(vtxundo);
    )

    bool WriteToDisk(CDiskBlockPos &pos, const uint256 &hashBlock)
    {
        // Open history file to append
        CAutoFile fileout = CAutoFile(OpenUndoFile(pos), SER_DISK, CLIENT_VERSION);
        if (!fileout)
            return error("CBlockUndo::WriteToDisk() : OpenUndoFile failed");

        // Write index header
        unsigned int nSize = fileout.GetSerializeSize(*this);
        fileout << FLATDATA(Params().MessageStart()) << nSize;

        // Write undo data
        long fileOutPos = ftell(fileout);
        if (fileOutPos < 0)
            return error("CBlockUndo::WriteToDisk() : ftell failed");
        pos.nPos = (unsigned int)fileOutPos;
        fileout << *this;

        // calculate & write checksum
        CHashWriter hasher(SER_GETHASH, PROTOCOL_VERSION);
        hasher << hashBlock;
        hasher << *this;
        fileout << hasher.GetHash();

        // Flush stdio buffers and commit to disk before returning
        fflush(fileout);
        if (!IsInitialBlockDownload())
            FileCommit(fileout);

        return true;
    }

    bool ReadFromDisk(const CDiskBlockPos &pos, const uint256 &hashBlock)
    {
        // Open history file to read
        CAutoFile filein = CAutoFile(OpenUndoFile(pos, true), SER_DISK, CLIENT_VERSION);
        if (!filein)
            return error("CBlockUndo::ReadFromDisk() : OpenBlockFile failed");

        // Read block
        uint256 hashChecksum;
        try {
            filein >> *this;
            filein >> hashChecksum;
        }
        catch (std::exception &e) {
            return error("%s() : deserialize or I/O error", __PRETTY_FUNCTION__);
        }

        // Verify checksum
        CHashWriter hasher(SER_GETHASH, PROTOCOL_VERSION);
        hasher << hashBlock;
        hasher << *this;
        if (hashChecksum != hasher.GetHash())
            return error("CBlockUndo::ReadFromDisk() : checksum mismatch");

        return true;
    }
};

/** pruned version of CTransaction: only retains metadata and unspent transaction outputs
 *
 * Serialized format:
 * - VARINT(nVersion)
 * - VARINT(nCode)
 * - unspentness bitvector, for vout[2] and further; least significant byte first
 * - the non-spent CTxOuts (via CTxOutCompressor)
 * - VARINT(nHeight)
 *
 * The nCode value consists of:
 * - bit 1: IsCoinBase()
 * - bit 2: vout[0] is not spent
 * - bit 4: vout[1] is not spent
 * - The higher bits encode N, the number of non-zero bytes in the following bitvector.
 *   - In case both bit 2 and bit 4 are unset, they encode N-1, as there must be at
 *     least one non-spent output).
 *
 * Example: 0104835800816115944e077fe7c803cfa57f29b36bf87c1d358bb85e
 *          <><><--------------------------------------------><---->
 *          |  \                  |                             /
 *    version   code             vout[1]                  height
 *
 *    - version = 1
 *    - code = 4 (vout[1] is not spent, and 0 non-zero bytes of bitvector follow)
 *    - unspentness bitvector: as 0 non-zero bytes follow, it has length 0
 *    - vout[1]: 835800816115944e077fe7c803cfa57f29b36bf87c1d35
 *               * 8358: compact amount representation for 60000000000 (600 BTC)
 *               * 00: special txout type pay-to-pubkey-hash
 *               * 816115944e077fe7c803cfa57f29b36bf87c1d35: address uint160
 *    - height = 203998
 *
 *
 * Example: 0109044086ef97d5790061b01caab50f1b8e9c50a5057eb43c2d9563a4eebbd123008c988f1a4a4de2161e0f50aac7f17e7f9555caa486af3b
 *          <><><--><--------------------------------------------------><----------------------------------------------><---->
 *         /  \   \                     |                                                           |                     /
 *  version  code  unspentness       vout[4]                                                     vout[16]           height
 *
 *  - version = 1
 *  - code = 9 (coinbase, neither vout[0] or vout[1] are unspent,
 *                2 (1, +1 because both bit 2 and bit 4 are unset) non-zero bitvector bytes follow)
 *  - unspentness bitvector: bits 2 (0x04) and 14 (0x4000) are set, so vout[2+2] and vout[14+2] are unspent
 *  - vout[4]: 86ef97d5790061b01caab50f1b8e9c50a5057eb43c2d9563a4ee
 *             * 86ef97d579: compact amount representation for 234925952 (2.35 BTC)
 *             * 00: special txout type pay-to-pubkey-hash
 *             * 61b01caab50f1b8e9c50a5057eb43c2d9563a4ee: address uint160
 *  - vout[16]: bbd123008c988f1a4a4de2161e0f50aac7f17e7f9555caa4
 *              * bbd123: compact amount representation for 110397 (0.001 BTC)
 *              * 00: special txout type pay-to-pubkey-hash
 *              * 8c988f1a4a4de2161e0f50aac7f17e7f9555caa4: address uint160
 *  - height = 120891
 */
class CCoins
{
public:
    // whether transaction is a coinbase
    bool fCoinBase;
    // whether transaction is a coinstake
    bool fCoinStake;

    // transaction timestamp
    int64_t nTime;

    // unspent transaction outputs; spent outputs are .IsNull(); spent outputs at the end of the array are dropped
    std::vector<CTxOut> vout;

    // at which height this transaction was included in the active blockchain
    int nHeight;

    // version of the CTransaction; accesses to this value should probably check for nHeight as well,
    // as new tx version will probably only be introduced at certain heights
    int nVersion;

// construct a CCoins from a CTransaction, at a given height
    CCoins(const CTransaction &tx, int nHeightIn) : fCoinBase(tx.IsCoinBase()), fCoinStake(tx.IsCoinStake()), nTime(tx.nTime), vout(tx.vout), nHeight(nHeightIn), nVersion(tx.nVersion) { }

    // empty constructor
    CCoins() : fCoinBase(false), fCoinStake(false), nTime(0), vout(0), nHeight(0), nVersion(0) { }


    // remove spent outputs at the end of vout
    void Cleanup() {
        while (vout.size() > 0 && vout.back().IsNull())
            vout.pop_back();
        if (vout.empty())
            std::vector<CTxOut>().swap(vout);
    }

    void swap(CCoins &to) {
        std::swap(to.fCoinBase, fCoinBase);
        std::swap(to.fCoinStake, fCoinStake);
        std::swap(to.nTime, nTime);
        to.vout.swap(vout);
        std::swap(to.nHeight, nHeight);
        std::swap(to.nVersion, nVersion);
    }

    // equality test
    friend bool operator==(const CCoins &a, const CCoins &b) {
        if (fDebug)
        {
            if (a.fCoinBase != b.fCoinBase)
                LogPrintf("CCoins: fCoinBase mismatch\n");
            if (a.fCoinStake != b.fCoinStake)
                LogPrintf("CCoins: fCoinStake mismatch\n");
            if (a.nTime != b.nTime)
                LogPrintf("CCoins: nTime mismatch\n");
            if (a.nHeight != b.nHeight)
                LogPrintf("CCoins: nHeight mismatch\n");
            if (a.nVersion != b.nVersion)
                LogPrintf("CCoins: nVersion mismatch\n");
            if (a.vout != b.vout)
                LogPrintf("CCoins: vout mismatch\n");
        }

        return a.fCoinBase == b.fCoinBase &&
                a.fCoinStake == b.fCoinStake &&
                a.nTime == b.nTime &&
                a.nHeight == b.nHeight &&
                a.nVersion == b.nVersion &&
                a.vout == b.vout;
    }
    friend bool operator!=(const CCoins &a, const CCoins &b) {
        return !(a == b);
    }

    // calculate number of bytes for the bitmask, and its number of non-zero bytes
    // each bit in the bitmask represents the availability of one output, but the
    // availabilities of the first two outputs are encoded separately
    void CalcMaskSize(unsigned int &nBytes, unsigned int &nNonzeroBytes) const {
        unsigned int nLastUsedByte = 0;
        for (unsigned int b = 0; 2+b*8 < vout.size(); b++) {
            bool fZero = true;
            for (unsigned int i = 0; i < 8 && 2+b*8+i < vout.size(); i++) {
                if (!vout[2+b*8+i].IsNull()) {
                    fZero = false;
                    continue;
                }
            }
            if (!fZero) {
                nLastUsedByte = b + 1;
                nNonzeroBytes++;
            }
        }
        nBytes += nLastUsedByte;
    }

    bool IsCoinBase() const {
        return fCoinBase;
    }
    
    bool IsCoinStake() const {
        return fCoinStake;
    }

    unsigned int GetSerializeSize(int nType, int nVersion) const {
        unsigned int nSize = 0;
        unsigned int nMaskSize = 0, nMaskCode = 0;
        CalcMaskSize(nMaskSize, nMaskCode);
        bool fFirst = vout.size() > 0 && !vout[0].IsNull();
        bool fSecond = vout.size() > 1 && !vout[1].IsNull();
        assert(fFirst || fSecond || nMaskCode);
        unsigned int nCode = 8*(nMaskCode - (fFirst || fSecond ? 0 : 1)) + (fCoinBase ? 1 : 0) + (fFirst ? 2 : 0) + (fSecond ? 4 : 0);
        // version
        nSize += ::GetSerializeSize(VARINT(this->nVersion), nType, nVersion);
        // size of header code
        nSize += ::GetSerializeSize(VARINT(nCode), nType, nVersion);
        // spentness bitmask
        nSize += nMaskSize;
        // txouts themself
        for (unsigned int i = 0; i < vout.size(); i++)
            if (!vout[i].IsNull())
                nSize += ::GetSerializeSize(CTxOutCompressor(REF(vout[i])), nType, nVersion);
        // height
        nSize += ::GetSerializeSize(VARINT(nHeight), nType, nVersion);
        // tx timestamp
        nSize += ::GetSerializeSize(VARINT(nTime), nType, nVersion);
        // coinstake
        nSize += 1;
        return nSize;
    }

    template<typename Stream>
    void Serialize(Stream &s, int nType, int nVersion) const {
        unsigned int nMaskSize = 0, nMaskCode = 0;
        CalcMaskSize(nMaskSize, nMaskCode);
        bool fFirst = vout.size() > 0 && !vout[0].IsNull();
        bool fSecond = vout.size() > 1 && !vout[1].IsNull();
        assert(fFirst || fSecond || nMaskCode);
        unsigned int nCode = 8*(nMaskCode - (fFirst || fSecond ? 0 : 1)) + (fCoinBase ? 1 : 0) + (fFirst ? 2 : 0) + (fSecond ? 4 : 0);
        // version
        ::Serialize(s, VARINT(this->nVersion), nType, nVersion);
        // header code
        ::Serialize(s, VARINT(nCode), nType, nVersion);
        // spentness bitmask
        for (unsigned int b = 0; b<nMaskSize; b++) {
            unsigned char chAvail = 0;
            for (unsigned int i = 0; i < 8 && 2+b*8+i < vout.size(); i++)
                if (!vout[2+b*8+i].IsNull())
                    chAvail |= (1 << i);
            ::Serialize(s, chAvail, nType, nVersion);
        }
        // txouts themself
        for (unsigned int i = 0; i < vout.size(); i++) {
            if (!vout[i].IsNull())
                ::Serialize(s, CTxOutCompressor(REF(vout[i])), nType, nVersion);
        }
        // coinbase height
        ::Serialize(s, VARINT(nHeight), nType, nVersion);
        // tx timestamp
        ::Serialize(s, VARINT(nTime), nType, nVersion);
        // coinstake
        unsigned char nCoinStake = fCoinStake ? 1 : 0;
        ::Serialize(s, nCoinStake, nType, nVersion);
    }

    template<typename Stream>
    void Unserialize(Stream &s, int nType, int nVersion) {
        unsigned int nCode = 0;
        // version
        ::Unserialize(s, VARINT(this->nVersion), nType, nVersion);
        // header code
        ::Unserialize(s, VARINT(nCode), nType, nVersion);
        fCoinBase = nCode & 1;
        std::vector<bool> vAvail(2, false);
        vAvail[0] = nCode & 2;
        vAvail[1] = nCode & 4;
        unsigned int nMaskCode = (nCode / 8) + ((nCode & 6) != 0 ? 0 : 1);
        // spentness bitmask
        while (nMaskCode > 0) {
            unsigned char chAvail = 0;
            ::Unserialize(s, chAvail, nType, nVersion);
            for (unsigned int p = 0; p < 8; p++) {
                bool f = (chAvail & (1 << p)) != 0;
                vAvail.push_back(f);
            }
            if (chAvail != 0)
                nMaskCode--;
        }
        // txouts themself
        vout.assign(vAvail.size(), CTxOut());
        for (unsigned int i = 0; i < vAvail.size(); i++) {
            if (vAvail[i])
                ::Unserialize(s, REF(CTxOutCompressor(vout[i])), nType, nVersion);
        }
        // coinbase height
        ::Unserialize(s, VARINT(nHeight), nType, nVersion);
        // tx timestamp
        ::Unserialize(s, VARINT(nTime), nType, nVersion);
        // coinstake
        unsigned char nCoinStake = 0;
        ::Unserialize(s, nCoinStake, nType, nVersion);
        fCoinStake = nCoinStake & 1;
        Cleanup();
    }

    // mark an outpoint spent, and construct undo information
    bool Spend(const COutPoint &out, CTxInUndo &undo) {
        if (out.n >= vout.size())
            return false;
        if (vout[out.n].IsNull())
            return false;
        undo = CTxInUndo(vout[out.n]);
        vout[out.n].SetNull();
        Cleanup();
        if (vout.size() == 0) {
            undo.nHeight = nHeight;
            undo.fCoinBase = fCoinBase;
            undo.nVersion = this->nVersion;
        }
        return true;
    }

    // mark a vout spent
    bool Spend(int nPos) {
        CTxInUndo undo;
        COutPoint out(0, nPos);
        return Spend(out, undo);
    }

    // check whether a particular output is still available
    bool IsAvailable(unsigned int nPos) const {
        return (nPos < vout.size() && !vout[nPos].IsNull());
    }

    // check whether the entire CCoins is spent
    // note that only !IsPruned() CCoins can be serialized
    bool IsPruned() const {
        BOOST_FOREACH(const CTxOut &out, vout)
            if (!out.IsNull())
                return false;
        return true;
    }
};

/** Closure representing one script verification
 *  Note that this stores references to the spending transaction */
class CScriptCheck
{
private:
    CScript scriptPubKey;
    const CTransaction *ptxTo;
    unsigned int nIn;
    unsigned int nFlags;
    unsigned int nBlockHeight;
    unsigned int nBlockTime;
    int nHashType;

public:
    CScriptCheck() {}
    CScriptCheck(const CCoins& txFromIn, const CTransaction& txToIn, unsigned int nInIn, unsigned int nFlagsIn, int nHashTypeIn, int nBlockHeight, int nBlockTime) :
        scriptPubKey(txFromIn.vout[txToIn.vin[nInIn].prevout.n].scriptPubKey),
        ptxTo(&txToIn), nIn(nInIn), nFlags(nFlagsIn), nHashType(nHashTypeIn) { }

    bool operator()() const;

    void swap(CScriptCheck &check) {
        scriptPubKey.swap(check.scriptPubKey);
        std::swap(ptxTo, check.ptxTo);
        std::swap(nIn, check.nIn);
        std::swap(nFlags, check.nFlags);
        std::swap(nHashType, check.nHashType);
        std::swap(nBlockHeight, check.nBlockHeight);
        std::swap(nBlockTime, check.nBlockTime);
    }
};



/** A transaction with a merkle branch linking it to the block chain. */
class CMerkleTx : public CTransaction
{
private:
    int GetDepthInMainChainINTERNAL(CBlockIndex* &pindexRet) const;
public:
    uint256 hashBlock;
    std::vector<uint256> vMerkleBranch;
    int nIndex;

    // memory only
    mutable bool fMerkleVerified;


    CMerkleTx()
    {
        Init();
    }

    CMerkleTx(const CTransaction& txIn) : CTransaction(txIn)
    {
        Init();
    }

    void Init()
    {
        hashBlock = 0;
        nIndex = -1;
        fMerkleVerified = false;
    }


    IMPLEMENT_SERIALIZE
    (
        nSerSize += SerReadWrite(s, *(CTransaction*)this, nType, nVersion, ser_action);
        nVersion = this->nVersion;
        READWRITE(hashBlock);
        READWRITE(vMerkleBranch);
        READWRITE(nIndex);
    )


    int SetMerkleBranch(const CBlock* pblock=NULL);

    // Return depth of transaction in blockchain:
    // -1  : not in blockchain, and not in memory pool (conflicted transaction)
    //  0  : in memory pool, waiting to be included in a block
    // >=1 : this many blocks deep in the main chain
    int GetDepthInMainChain(CBlockIndex* &pindexRet) const;
    int GetDepthInMainChain() const { CBlockIndex *pindexRet; return GetDepthInMainChain(pindexRet); }
    bool IsInMainChain() const { CBlockIndex *pindexRet; return GetDepthInMainChainINTERNAL(pindexRet) > 0; }
    int GetBlocksToMaturity() const;
    bool AcceptToMemoryPool(bool fCheckInputs=true, bool fLimitFree=true);
};



/** Data structure that represents a partial merkle tree.
 *
 * It respresents a subset of the txid's of a known block, in a way that
 * allows recovery of the list of txid's and the merkle root, in an
 * authenticated way.
 *
 * The encoding works as follows: we traverse the tree in depth-first order,
 * storing a bit for each traversed node, signifying whether the node is the
 * parent of at least one matched leaf txid (or a matched txid itself). In
 * case we are at the leaf level, or this bit is 0, its merkle node hash is
 * stored, and its children are not explorer further. Otherwise, no hash is
 * stored, but we recurse into both (or the only) child branch. During
 * decoding, the same depth-first traversal is performed, consuming bits and
 * hashes as they written during encoding.
 *
 * The serialization is fixed and provides a hard guarantee about the
 * encoded size:
 *
 *   SIZE <= 10 + ceil(32.25*N)
 *
 * Where N represents the number of leaf nodes of the partial tree. N itself
 * is bounded by:
 *
 *   N <= total_transactions
 *   N <= 1 + matched_transactions*tree_height
 *
 * The serialization format:
 *  - uint32     total_transactions (4 bytes)
 *  - varint     number of hashes   (1-3 bytes)
 *  - uint256[]  hashes in depth-first order (<= 32*N bytes)
 *  - varint     number of bytes of flag bits (1-3 bytes)
 *  - byte[]     flag bits, packed per 8 in a byte, least significant bit first (<= 2*N-1 bits)
 * The size constraints follow from this.
 */
class CPartialMerkleTree
{
protected:
    // the total number of transactions in the block
    unsigned int nTransactions;

    // node-is-parent-of-matched-txid bits
    std::vector<bool> vBits;

    // txids and internal hashes
    std::vector<uint256> vHash;

    // flag set when encountering invalid data
    bool fBad;

    // helper function to efficiently calculate the number of nodes at given height in the merkle tree
    unsigned int CalcTreeWidth(int height) {
        return (nTransactions+(1 << height)-1) >> height;
    }

    // calculate the hash of a node in the merkle tree (at leaf level: the txid's themself)
    uint256 CalcHash(int height, unsigned int pos, const std::vector<uint256> &vTxid);

    // recursive function that traverses tree nodes, storing the data as bits and hashes
    void TraverseAndBuild(int height, unsigned int pos, const std::vector<uint256> &vTxid, const std::vector<bool> &vMatch);

    // recursive function that traverses tree nodes, consuming the bits and hashes produced by TraverseAndBuild.
    // it returns the hash of the respective node.
    uint256 TraverseAndExtract(int height, unsigned int pos, unsigned int &nBitsUsed, unsigned int &nHashUsed, std::vector<uint256> &vMatch);

public:

    // serialization implementation
    IMPLEMENT_SERIALIZE(
        READWRITE(nTransactions);
        READWRITE(vHash);
        std::vector<unsigned char> vBytes;
        if (fRead) {
            READWRITE(vBytes);
            CPartialMerkleTree &us = *(const_cast<CPartialMerkleTree*>(this));
            us.vBits.resize(vBytes.size() * 8);
            for (unsigned int p = 0; p < us.vBits.size(); p++)
                us.vBits[p] = (vBytes[p / 8] & (1 << (p % 8))) != 0;
            us.fBad = false;
        } else {
            vBytes.resize((vBits.size()+7)/8);
            for (unsigned int p = 0; p < vBits.size(); p++)
                vBytes[p / 8] |= vBits[p] << (p % 8);
            READWRITE(vBytes);
        }
    )

    // Construct a partial merkle tree from a list of transaction id's, and a mask that selects a subset of them
    CPartialMerkleTree(const std::vector<uint256> &vTxid, const std::vector<bool> &vMatch);

    CPartialMerkleTree();

    // extract the matching txid's represented by this partial merkle tree.
    // returns the merkle root, or 0 in case of failure
    uint256 ExtractMatches(std::vector<uint256> &vMatch);
};


/** Nodes collect new transactions into a block, hash them into a hash tree,
 * and scan through nonce values to make the block's hash satisfy proof-of-work
 * requirements.  When they solve the proof-of-work, they broadcast the block
 * to everyone and the block is added to the block chain.  The first transaction
 * in the block is a special one that creates a new coin owned by the creator
 * of the block.
 */
class CBlockHeader
{
public:
    // header
    static const int CURRENT_VERSION=7;
    int nVersion;
    uint256 hashPrevBlock;
    uint256 hashMerkleRoot;
    unsigned int nTime;
    unsigned int nBits;
    unsigned int nNonce;

    CBlockHeader()
    {
        SetNull();
    }

    IMPLEMENT_SERIALIZE
    (
        READWRITE(this->nVersion);
        nVersion = this->nVersion;
        READWRITE(hashPrevBlock);
        READWRITE(hashMerkleRoot);
        READWRITE(nTime);
        READWRITE(nBits);
        READWRITE(nNonce);
    )

    void SetNull()
    {
        nVersion = CBlockHeader::CURRENT_VERSION;
        hashPrevBlock = 0;
        hashMerkleRoot = 0;
        nTime = 0;
        nBits = 0;
        nNonce = 0;
    }

    bool IsNull() const
    {
        return (nBits == 0);
    }

    uint256 GetHash() const
    {
        if (nVersion > 6)
            return Hash(BEGIN(nVersion), END(nNonce));
        else
            return GetPoWHash();
    }

    uint256 GetPoWHash() const
    {
        return scrypt_blockhash(CVOIDBEGIN(nVersion));
    }

    int64_t GetBlockTime() const
    {
        return (int64)nTime;
    }

    void UpdateTime(const CBlockIndex* pindexPrev);
};


class CBlock : public CBlockHeader
{
public:
    // network and disk
    std::vector<CTransaction> vtx;

    // ppcoin: block signature - signed by one of the coin base txout[N]'s owner
    std::vector<unsigned char> vchBlockSig;

    // memory only
    mutable std::vector<uint256> vMerkleTree;

    CBlock()
    {
        SetNull();
    }

    CBlock(const CBlockHeader &header)
    {
        SetNull();
        *((CBlockHeader*)this) = header;
    }

    IMPLEMENT_SERIALIZE
    (
        READWRITE(*(CBlockHeader*)this);
            // ConnectBlock depends on vtx following header to generate CDiskTxPos
        if (!(nType & (SER_GETHASH|SER_BLOCKHEADERONLY)))
        {
            READWRITE(vtx);
            READWRITE(vchBlockSig);
        }
        else if (fRead)
        {
            const_cast<CBlock*>(this)->vtx.clear();
            const_cast<CBlock*>(this)->vchBlockSig.clear();
        }
    )

    void SetNull()
    {
        CBlockHeader::SetNull();
        vtx.clear();
        vMerkleTree.clear();
        vchBlockSig.clear(); // ppcoin
    }

    uint256 GetPoWHash() const
    {
        return scrypt_blockhash(CVOIDBEGIN(nVersion));
    }

    // ppcoin: entropy bit for stake modifier if chosen by modifier
    unsigned int GetStakeEntropyBit() const
    {
        // Take last bit of block hash as entropy bit
        unsigned int nEntropyBit = ((GetHash().Get64()) & 1llu);
        if (fDebug && GetBoolArg("-printstakemodifier", false))
            LogPrintf("GetStakeEntropyBit: hashBlock=%s nEntropyBit=%u\n", GetHash().ToString().c_str(), nEntropyBit);
        return nEntropyBit;
    }

    // ppcoin: two types of block: proof-of-work or proof-of-stake
    bool IsProofOfStake() const
    {
        return (vtx.size() > 1 && vtx[1].IsCoinStake());
    }

    bool IsProofOfWork() const
    {
        return !IsProofOfStake();
    }

    std::pair<COutPoint, unsigned int> GetProofOfStake() const
    {
        return IsProofOfStake()? std::make_pair(vtx[1].vin[0].prevout, vtx[1].nTime) : std::make_pair(COutPoint(), (unsigned int)0);
    }

    // ppcoin: get max transaction timestamp
    int64_t GetMaxTransactionTime() const
    {
        int64_t maxTransactionTime = 0;
        BOOST_FOREACH(const CTransaction& tx, vtx)
            maxTransactionTime = std::max(maxTransactionTime, (int64_t)tx.nTime);
        return maxTransactionTime;
    }

    CBlockHeader GetBlockHeader() const
    {
        CBlockHeader block;
        block.nVersion       = nVersion;
        block.hashPrevBlock  = hashPrevBlock;
        block.hashMerkleRoot = hashMerkleRoot;
        block.nTime          = nTime;
        block.nBits          = nBits;
        block.nNonce         = nNonce;
        return block;
    }

    uint256 BuildMerkleTree() const
    {
        vMerkleTree.clear();
        BOOST_FOREACH(const CTransaction& tx, vtx)
            vMerkleTree.push_back(tx.GetHash());
        int j = 0;
        for (int nSize = vtx.size(); nSize > 1; nSize = (nSize + 1) / 2)
        {
            for (int i = 0; i < nSize; i += 2)
            {
                int i2 = std::min(i+1, nSize-1);
                vMerkleTree.push_back(Hash(BEGIN(vMerkleTree[j+i]),  END(vMerkleTree[j+i]),
                                           BEGIN(vMerkleTree[j+i2]), END(vMerkleTree[j+i2])));
            }
            j += nSize;
        }
        return (vMerkleTree.empty() ? 0 : vMerkleTree.back());
    }

    const uint256 &GetTxHash(unsigned int nIndex) const {
        assert(vMerkleTree.size() > 0); // BuildMerkleTree must have been called first
        assert(nIndex < vtx.size());
        return vMerkleTree[nIndex];
    }

    std::vector<uint256> GetMerkleBranch(int nIndex) const
    {
        if (vMerkleTree.empty())
            BuildMerkleTree();
        std::vector<uint256> vMerkleBranch;
        int j = 0;
        for (int nSize = vtx.size(); nSize > 1; nSize = (nSize + 1) / 2)
        {
            int i = std::min(nIndex^1, nSize-1);
            vMerkleBranch.push_back(vMerkleTree[j+i]);
            nIndex >>= 1;
            j += nSize;
        }
        return vMerkleBranch;
    }

    static uint256 CheckMerkleBranch(uint256 hash, const std::vector<uint256>& vMerkleBranch, int nIndex)
    {
        if (nIndex == -1)
            return 0;
        BOOST_FOREACH(const uint256& otherside, vMerkleBranch)
        {
            if (nIndex & 1)
                hash = Hash(BEGIN(otherside), END(otherside), BEGIN(hash), END(hash));
            else
                hash = Hash(BEGIN(hash), END(hash), BEGIN(otherside), END(otherside));
            nIndex >>= 1;
        }
        return hash;
    }

    bool WriteToDisk(CDiskBlockPos &pos)
    {
        // Open history file to append
        CAutoFile fileout = CAutoFile(OpenBlockFile(pos), SER_DISK, CLIENT_VERSION);
        if (!fileout)
            return error("CBlock::WriteToDisk() : OpenBlockFile failed");

        // Write index header
        unsigned int nSize = fileout.GetSerializeSize(*this);
        fileout << FLATDATA(Params().MessageStart()) << nSize;

        // Write block
        long fileOutPos = ftell(fileout);
        if (fileOutPos < 0)
            return error("CBlock::WriteToDisk() : ftell failed");
        pos.nPos = (unsigned int)fileOutPos;
        fileout << *this;

        // Flush stdio buffers and commit to disk before returning
        fflush(fileout);
        if (!IsInitialBlockDownload())
            FileCommit(fileout);

        return true;
    }

    bool ReadFromDisk(const CDiskBlockPos &pos)
    {
        SetNull();

        // Open history file to read
        CAutoFile filein = CAutoFile(OpenBlockFile(pos, true), SER_DISK, CLIENT_VERSION);
        if (!filein)
            return error("CBlock::ReadFromDisk() : OpenBlockFile failed");

        // Read block
        try {
            filein >> *this;
        }
        catch (std::exception &e) {
            return error("%s() : deserialize or I/O error", __PRETTY_FUNCTION__);
        }

        // Check the header
        if (IsProofOfWork() && !CheckProofOfWork(GetPoWHash(), nBits))
            return error("CBlock::ReadFromDisk() : errors in block header");

        return true;
    }

    void print() const
    {
        LogPrintf("CBlock(hash=%s, input=%s, PoW=%s, ver=%d, hashPrevBlock=%s, hashMerkleRoot=%s, nTime=%u, nBits=%08x, nNonce=%u, vtx=%"PRIszu", vchBlockSig=%s)\n",
            GetHash().ToString().c_str(),
            HexStr(BEGIN(nVersion),BEGIN(nVersion)+80,false).c_str(),
            GetPoWHash().ToString().c_str(),
            nVersion,
            hashPrevBlock.ToString().c_str(),
            hashMerkleRoot.ToString().c_str(),
            nTime, nBits, nNonce,
            vtx.size(),
            HexStr(vchBlockSig.begin(), vchBlockSig.end()).c_str());
        for (unsigned int i = 0; i < vtx.size(); i++)
        {
            LogPrintf("  ");
            vtx[i].print();
        }
        LogPrintf("  vMerkleTree: ");
        for (unsigned int i = 0; i < vMerkleTree.size(); i++)
            LogPrintf("%s ", vMerkleTree[i].ToString().c_str());
        LogPrintf("\n");
    }


    /** Undo the effects of this block (with given index) on the UTXO set represented by coins.
     *  In case pfClean is provided, operation will try to be tolerant about errors, and *pfClean
     *  will be true if no problems were found. Otherwise, the return value will be false in case
     *  of problems. Note that in any case, coins may be modified. */
    bool DisconnectBlock(CValidationState &state, CBlockIndex *pindex, CCoinsViewCache &coins, bool *pfClean = NULL);

    // Apply the effects of this block (with given index) on the UTXO set represented by coins
    bool ConnectBlock(CValidationState &state, CBlockIndex *pindex, CCoinsViewCache &coins, bool fJustCheck=false);

    // Read a block from disk
    bool ReadFromDisk(const CBlockIndex* pindex);

    // Add this block to the block index, and if necessary, switch the active block chain to this
    bool AddToBlockIndex(CValidationState &state, const CDiskBlockPos &pos, const uint256 &hashProof);

    // Context-independent validity checks
    bool CheckBlock(CValidationState &state, bool fCheckPOW=true, bool fCheckMerkleRoot=true, bool fCheckSig=true, int nHeight=0) const;

    // Store block on disk
    // if dbp is provided, the file is known to already reside on disk
    bool AcceptBlock(CValidationState &state, CDiskBlockPos *dbp = NULL);

    // ppcoin: calculate total coin age spent in block
    bool GetCoinAge(uint64_t& nCoinAge) const;
    bool SignBlock(CWallet& keystore, int64_t nFees);
    bool CheckBlockSignature() const;
};

enum BlockStatus {
    BLOCK_VALID_UNKNOWN      =    0,
    BLOCK_VALID_HEADER       =    1, // parsed, version ok, hash satisfies claimed PoW, 1 <= vtx count <= max, timestamp not in future
    BLOCK_VALID_TREE         =    2, // parent found, difficulty matches, timestamp >= median previous, checkpoint
    BLOCK_VALID_TRANSACTIONS =    3, // only first tx is coinbase, 2 <= coinbase input script length <= 100, transactions valid, no duplicate txids, sigops, size, merkle root
    BLOCK_VALID_CHAIN        =    4, // outputs do not overspend inputs, no double spends, coinbase output ok, immature coinbase spends, BIP30
    BLOCK_VALID_SCRIPTS      =    5, // scripts/signatures ok
    BLOCK_VALID_MASK         =    7,

    BLOCK_HAVE_DATA          =    8, // full block available in blk*.dat
    BLOCK_HAVE_UNDO          =   16, // undo data available in rev*.dat
    BLOCK_HAVE_MASK          =   24,

    BLOCK_FAILED_VALID       =   32, // stage after last reached validness failed
    BLOCK_FAILED_CHILD       =   64, // descends from failed block
    BLOCK_FAILED_MASK        =   96
};


/** The block chain is a tree shaped structure starting with the
 * genesis block at the root, with each block potentially having multiple
 * candidates to be the next block.  pprev and pnext link a path through the
 * main/longest chain.  A blockindex may have multiple pprev pointing back
 * to it, but pnext will only point forward to the longest branch, or will
 * be null if the block is not part of the longest chain.
 */
class CClamour
{
public:
    int nHeight;
    uint256 txid;
    std::string strHash;
    std::string strURL;

    CClamour()
    {
    }

    CClamour(int nHeightIn, uint256 txidIn, std::string& strHashIn, std::string& strURLIn)
    {
        nHeight = nHeightIn;
        txid = txidIn;
        strHash = strHashIn;
        strURL = strURLIn;
    }

    IMPLEMENT_SERIALIZE
    (
        READWRITE(nHeight);
        READWRITE(txid);
        READWRITE(strHash);
        READWRITE(strURL);
    )
};

class CBlockIndex
{
public:
    const uint256* phashBlock;
    CBlockIndex* pprev;
    CBlockIndex* pnext;
    unsigned int nFile;
    // Byte offset within blk?????.dat where this block's data is stored
    unsigned int nDataPos;

    // Byte offset within rev?????.dat where this block's undo data is stored
    unsigned int nUndoPos;
    // (memory only) Total amount of work (expected number of hashes) in the chain up to and including this block
    uint256 nChainWork;
    int nHeight;
    // Number of transactions in this block.
    // Note: in a potential headers-first mode, this number cannot be relied upon
    unsigned int nTx;
    // (memory only) Number of transactions in the chain up to and including this block
    unsigned int nChainTx; // change to 64-bit type when necessary; won't happen before 2030

     // Verification status of this block. See enum BlockStatus
    unsigned int nStatus;

    int64_t nMint;
    int64_t nMoneySupply;
    int64_t nDigsupply;
    int64_t nStakeSupply;

    std::vector<CClamour> vClamour;

    unsigned int nFlags;  // ppcoin: block index flags
    enum  
    {
        BLOCK_PROOF_OF_STAKE = (1 << 0), // is proof-of-stake block
        BLOCK_STAKE_ENTROPY  = (1 << 1), // entropy bit for stake modifier
        BLOCK_STAKE_MODIFIER = (1 << 2), // regenerated stake modifier
    };

    uint64_t nStakeModifier; // hash modifier for proof-of-stake

    // proof-of-stake specific fields
    COutPoint prevoutStake;
    unsigned int nStakeTime;

    uint256 hashProof;

    // block header
    int nVersion;
    uint256 hashMerkleRoot;
    unsigned int nTime;
    unsigned int nBits;
    unsigned int nNonce;

    mutable bool fSupportChecked; // did we check the speech of the staking transaction for 'clamour' support yet?
    mutable std::set<std::string> setSupport; // CLAMour pid strings supported by this block

    CBlockIndex()
    {
        phashBlock = NULL;
        pprev = NULL;
        pnext = NULL;
        nFile = 0;
        nDataPos = 0;
        nUndoPos = 0;
        nHeight = 0;
        nChainWork = 0;
        nHeight = 0;
        nMint = 0;
        nMoneySupply = 0;
        nDigsupply = 0;
        nStakeSupply = 0;
        nFlags = 0;
        nStakeModifier = 0;
        hashProof = 0;
        prevoutStake.SetNull();
        nStakeTime = 0;
        nTx = 0;
        nChainTx = 0;
        nStatus = 0;

        nVersion       = 0;
        hashMerkleRoot = 0;
        nTime          = 0;
        nBits          = 0;
        nNonce         = 0;

        fSupportChecked = false;
    }

    CBlockIndex(CBlock& block)
    {
        phashBlock = NULL;
        pprev = NULL;
        pnext = NULL;
        nFile = 0;
        nDataPos = 0;
        nUndoPos = 0;
        nTx = 0;
        nChainTx = 0;
        nStatus = 0;
        nHeight = 0;
        nChainWork = 0;
        nMint = 0;
        nMoneySupply = 0;
        nDigsupply = 0;
        nStakeSupply = 0;
        nFlags = 0;
        nStakeModifier = 0;
        hashProof = 0;
        if (block.IsProofOfStake())
        {
            SetProofOfStake();
            prevoutStake = block.vtx[1].vin[0].prevout;
            nStakeTime = block.vtx[1].nTime;
        }
        else
        {
            prevoutStake.SetNull();
            nStakeTime = 0;
        }

        nVersion       = block.nVersion;
        hashMerkleRoot = block.hashMerkleRoot;
        nTime          = block.nTime;
        nBits          = block.nBits;
        nNonce         = block.nNonce;

        fSupportChecked = false;
    }

    CDiskBlockPos GetBlockPos() const {
        CDiskBlockPos ret;
        if (nStatus & BLOCK_HAVE_DATA) {
            ret.nFile = nFile;
            ret.nPos  = nDataPos;
        } else
            ret.SetNull();
        return ret;
    }

    CDiskBlockPos GetUndoPos() const {
        CDiskBlockPos ret;
        if (nStatus & BLOCK_HAVE_UNDO) {
            ret.nFile = nFile;
            ret.nPos  = nUndoPos;
        } else
            ret.SetNull();
        return ret;
    }

    CBlock GetBlockHeader() const
    {
        CBlock block;
        block.nVersion       = nVersion;
        if (pprev)
            block.hashPrevBlock = pprev->GetBlockHash();
        block.hashMerkleRoot = hashMerkleRoot;
        block.nTime          = nTime;
        block.nBits          = nBits;
        block.nNonce         = nNonce;
        return block;
    }

    uint256 GetBlockHash() const
    {
        return *phashBlock;
    }

    int64_t GetBlockTime() const
    {
        return (int64_t)nTime;
    }

    CBigNum GetBlockWork() const
    {
        CBigNum bnTarget;
        bnTarget.SetCompact(nBits);
        if (bnTarget <= 0)
            return 0;
        return (CBigNum(1)<<256) / (bnTarget+1);
    }


    bool IsInMainChain() const
    {
        return (pnext || this == pindexBest);
    }

    bool CheckIndex() const
    {
        return true;
    }

    int64_t GetPastTimeLimit() const
    {
        if (IsProtocolV2(nHeight))
            return GetBlockTime();
        else
            return GetMedianTimePast();
    }

    enum { nMedianTimeSpan=11 };

    int64_t GetMedianTimePast() const
    {
        int64_t pmedian[nMedianTimeSpan];
        int64_t* pbegin = &pmedian[nMedianTimeSpan];
        int64_t* pend = &pmedian[nMedianTimeSpan];

        const CBlockIndex* pindex = this;
        for (int i = 0; i < nMedianTimeSpan && pindex; i++, pindex = pindex->pprev)
            *(--pbegin) = pindex->GetBlockTime();

        std::sort(pbegin, pend);
        return pbegin[(pend - pbegin)/2];
    }

    /**
     * Returns true if there are nRequired or more blocks of minVersion or above
     * in the last nToCheck blocks, starting at pstart and going backwards.
     */
    static bool IsSuperMajority(int minVersion, const CBlockIndex* pstart,
                                unsigned int nRequired, unsigned int nToCheck);


    bool IsProofOfWork() const
    {
        return !(nFlags & BLOCK_PROOF_OF_STAKE);
    }

    bool IsProofOfStake() const
    {
        return (nFlags & BLOCK_PROOF_OF_STAKE);
    }

    void SetProofOfStake()
    {
        nFlags |= BLOCK_PROOF_OF_STAKE;
    }

    unsigned int GetStakeEntropyBit() const
    {
        return ((nFlags & BLOCK_STAKE_ENTROPY) >> 1);
    }

    bool SetStakeEntropyBit(unsigned int nEntropyBit)
    {
        if (nEntropyBit > 1)
            return false;
        nFlags |= (nEntropyBit? BLOCK_STAKE_ENTROPY : 0);
        return true;
    }

    bool GeneratedStakeModifier() const
    {
        return (nFlags & BLOCK_STAKE_MODIFIER);
    }

    void SetStakeModifier(uint64_t nModifier, bool fGeneratedStakeModifier)
    {
        nStakeModifier = nModifier;
        if (fGeneratedStakeModifier)
            nFlags |= BLOCK_STAKE_MODIFIER;
    }

    std::set<std::string> GetSupport() const;

    std::string ToString() const
    {
        return strprintf("CBlockIndex(pprev=%p, pnext=%p, nHeight=%d, nMint=%s, nMoneySupply=%s, nFlags=(%s)(%d)(%s), nStakeModifier=%016, hashProof=%s, prevoutStake=(%s), nStakeTime=%d merkle=%s, hashBlock=%s)",
            pprev, pnext, nHeight,
            FormatMoney(nMint).c_str(), FormatMoney(nMoneySupply).c_str(),
            GeneratedStakeModifier() ? "MOD" : "-", GetStakeEntropyBit(), IsProofOfStake()? "PoS" : "PoW",
            nStakeModifier,
            hashProof.ToString().c_str(),
            prevoutStake.ToString().c_str(), nStakeTime,
            hashMerkleRoot.ToString().c_str(),
            GetBlockHash().ToString().c_str());
    }

    void print() const
    {
        LogPrintf("%s\n", ToString().c_str());
    }
};

struct CBlockIndexWorkComparator
{
    bool operator()(CBlockIndex *pa, CBlockIndex *pb) {
        if (pa->nChainWork > pb->nChainWork) return false;
        if (pa->nChainWork < pb->nChainWork) return true;

        if (pa->GetBlockHash() < pb->GetBlockHash()) return false;
        if (pa->GetBlockHash() > pb->GetBlockHash()) return true;

        return false; // identical blocks
    }
};

/** Used to marshal pointers into hashes for db storage. */
class CDiskBlockIndex : public CBlockIndex
{
private:
    uint256 blockHash;

public:
    uint256 hashPrev;
    uint256 hashNext;

    CDiskBlockIndex()
    {
        hashPrev = 0;
        hashNext = 0;
        blockHash = 0;
    }

    explicit CDiskBlockIndex(CBlockIndex* pindex) : CBlockIndex(*pindex)
    {
        hashPrev = (pprev ? pprev->GetBlockHash() : 0);
        hashNext = (pnext ? pnext->GetBlockHash() : 0);
    }

    IMPLEMENT_SERIALIZE
    (
        if (!(nType & SER_GETHASH))
            READWRITE(nVersion);

        READWRITE(nHeight);
        READWRITE(VARINT(nStatus));
        READWRITE(VARINT(nTx));
        if (nStatus & (BLOCK_HAVE_DATA | BLOCK_HAVE_UNDO))
            READWRITE(VARINT(nFile));
        if (nStatus & BLOCK_HAVE_DATA)
            READWRITE(VARINT(nDataPos));
        if (nStatus & BLOCK_HAVE_UNDO)
            READWRITE(VARINT(nUndoPos));

        // block header
        READWRITE(this->nVersion);
        READWRITE(hashPrev);
        READWRITE(hashMerkleRoot);
        READWRITE(nTime);
        READWRITE(nBits);
        READWRITE(nNonce);
        READWRITE(blockHash);

        //clam
        READWRITE(hashNext);
        READWRITE(nMint);
        READWRITE(nMoneySupply);
        READWRITE(nDigsupply);
        READWRITE(nStakeSupply);
        READWRITE(vClamour);
        READWRITE(nFlags);
        READWRITE(nStakeModifier);

        if (IsProofOfStake())
        {
            READWRITE(prevoutStake);
            READWRITE(nStakeTime);
        }
        else if (fRead)
        {
            const_cast<CDiskBlockIndex*>(this)->prevoutStake.SetNull();
            const_cast<CDiskBlockIndex*>(this)->nStakeTime = 0;
        }
        READWRITE(hashProof);


    )

    uint256 GetBlockHash() const
    {
        if (fUseFastIndex && (nTime < GetAdjustedTime() - 24 * 60 * 60) && blockHash != 0)
            return blockHash;

        CBlock block;
        block.nVersion        = nVersion;
        block.hashPrevBlock   = hashPrev;
        block.hashMerkleRoot  = hashMerkleRoot;
        block.nTime           = nTime;
        block.nBits           = nBits;
        block.nNonce          = nNonce;

        const_cast<CDiskBlockIndex*>(this)->blockHash = block.GetHash();

        return blockHash;
    }

    std::string ToString() const
    {
        std::string str = "CDiskBlockIndex(";
        str += CBlockIndex::ToString();
        str += strprintf("\n                hashBlock=%s, hashPrev=%s, hashNext=%s)",
            GetBlockHash().ToString(),
            hashPrev.ToString(),
            hashNext.ToString());
        return str;
    }
};



/** Capture information about block/transaction validation */
class CValidationState {
private:
    enum mode_state {
        MODE_VALID,   // everything ok
        MODE_INVALID, // network rule violation (DoS value may be set)
        MODE_ERROR,   // run-time error
    } mode;
    int nDoS;
    bool corruptionPossible;
public:
    CValidationState() : mode(MODE_VALID), nDoS(0), corruptionPossible(false) {}
    bool DoS(int level, bool ret = false, bool corruptionIn = false) {
        if (mode == MODE_ERROR)
            return ret;
        nDoS += level;
        mode = MODE_INVALID;
        corruptionPossible = corruptionIn;
        return ret;
    }
    bool Invalid(bool ret = false) {
        return DoS(0, ret);
    }
    bool Error() {
        mode = MODE_ERROR;
        return false;
    }
    bool Abort(const std::string &msg) {
        AbortNode(msg);
        return Error();
    }
    bool IsValid() {
        return mode == MODE_VALID;
    }
    bool IsInvalid() {
        return mode == MODE_INVALID;
    }
    bool IsError() {
        return mode == MODE_ERROR;
    }
    bool IsInvalid(int &nDoSOut) {
        if (IsInvalid()) {
            nDoSOut = nDoS;
            return true;
        }
        return false;
    }
    bool CorruptionPossible() {
        return corruptionPossible;
    }
};

class CTxMemPool
{
public:
    mutable CCriticalSection cs;
    std::map<uint256, CTransaction> mapTx;
    std::map<COutPoint, CInPoint> mapNextTx;

    bool accept(CValidationState &state, CTransaction &tx, bool fCheckInputs, bool fLimitFree, bool* pfMissingInputs, bool fRejectInsaneFee = false);
    bool addUnchecked(const uint256& hash, const CTransaction &tx);
    bool remove(const CTransaction &tx, bool fRecursive = false);
    bool removeConflicts(const CTransaction &tx);
    void clear();
    void queryHashes(std::vector<uint256>& vtxid);
    void pruneSpent(const uint256& hash, CCoins &coins);

    unsigned long size()
    {
        LOCK(cs);
        return mapTx.size();
    }

    bool exists(uint256 hash)
    {
        return (mapTx.count(hash) != 0);
    }

    CTransaction& lookup(uint256 hash)
    {
        return mapTx[hash];
    }
};

extern CTxMemPool mempool;


/** Describes a place in the block chain to another node such that if the
 * other node doesn't have the same branch, it can find a recent common trunk.
 * The further back it is, the further before the fork it may be.
 */
class CBlockLocator
{
protected:
    std::vector<uint256> vHave;
public:

    CBlockLocator()
    {
    }

    explicit CBlockLocator(const CBlockIndex* pindex)
    {
        Set(pindex);
    }

    explicit CBlockLocator(uint256 hashBlock)
    {
        std::map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashBlock);
        if (mi != mapBlockIndex.end())
            Set((*mi).second);
    }

    CBlockLocator(const std::vector<uint256>& vHaveIn)
    {
        vHave = vHaveIn;
    }

    IMPLEMENT_SERIALIZE
    (
        if (!(nType & SER_GETHASH))
            READWRITE(nVersion);
        READWRITE(vHave);
    )

    void SetNull()
    {
        vHave.clear();
    }

    bool IsNull()
    {
        return vHave.empty();
    }

    void Set(const CBlockIndex* pindex)
    {
        vHave.clear();
        int nStep = 1;
        while (pindex)
        {
            vHave.push_back(pindex->GetBlockHash());

            // Exponentially larger steps back
            for (int i = 0; pindex && i < nStep; i++)
                pindex = pindex->pprev;
            if (vHave.size() > 10)
                nStep *= 2;
        }
        vHave.push_back(Params().HashGenesisBlock());
    }

    int GetDistanceBack()
    {
        // Retrace how far back it was in the sender's branch
        int nDistance = 0;
        int nStep = 1;
        BOOST_FOREACH(const uint256& hash, vHave)
        {
            std::map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hash);
            if (mi != mapBlockIndex.end())
            {
                CBlockIndex* pindex = (*mi).second;
                if (pindex->IsInMainChain())
                    return nDistance;
            }
            nDistance += nStep;
            if (nDistance > 10)
                nStep *= 2;
        }
        return nDistance;
    }

    CBlockIndex* GetBlockIndex()
    {
        // Find the first block the caller has in the main chain
        BOOST_FOREACH(const uint256& hash, vHave)
        {
            std::map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hash);
            if (mi != mapBlockIndex.end())
            {
                CBlockIndex* pindex = (*mi).second;
                if (pindex->IsInMainChain())
                    return pindex;
            }
        }
        return pindexGenesisBlock;
    }

    uint256 GetBlockHash()
    {
        // Find the first block the caller has in the main chain
        BOOST_FOREACH(const uint256& hash, vHave)
        {
            std::map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hash);
            if (mi != mapBlockIndex.end())
            {
                CBlockIndex* pindex = (*mi).second;
                if (pindex->IsInMainChain())
                    return hash;
            }
        }
        return Params().HashGenesisBlock();
    }

    int GetHeight()
    {
        CBlockIndex* pindex = GetBlockIndex();
        if (!pindex)
            return 0;
        return pindex->nHeight;
    }
};




class CBlockFileInfo
{
public:
    unsigned int nBlocks;      // number of blocks stored in file
    unsigned int nSize;        // number of used bytes of block file
    unsigned int nUndoSize;    // number of used bytes in the undo file
    unsigned int nHeightFirst; // lowest height of block in file
    unsigned int nHeightLast;  // highest height of block in file
    uint64_t nTimeFirst;         // earliest time of block in file
    uint64_t nTimeLast;          // latest time of block in file

    IMPLEMENT_SERIALIZE(
        READWRITE(VARINT(nBlocks));
        READWRITE(VARINT(nSize));
        READWRITE(VARINT(nUndoSize));
        READWRITE(VARINT(nHeightFirst));
        READWRITE(VARINT(nHeightLast));
        READWRITE(VARINT(nTimeFirst));
        READWRITE(VARINT(nTimeLast));
     )

     void SetNull() {
         nBlocks = 0;
         nSize = 0;
         nUndoSize = 0;
         nHeightFirst = 0;
         nHeightLast = 0;
         nTimeFirst = 0;
         nTimeLast = 0;
     }

     CBlockFileInfo() {
         SetNull();
     }

     std::string ToString() const {
         return strprintf("CBlockFileInfo(blocks=%u, size=%u, heights=%u..%u, time=%s..%s)", nBlocks, nSize, nHeightFirst, nHeightLast, DateTimeStrFormat("%Y-%m-%d", nTimeFirst).c_str(), DateTimeStrFormat("%Y-%m-%d", nTimeLast).c_str());
     }

     // update statistics (does not update nSize)
     void AddBlock(unsigned int nHeightIn, uint64_t nTimeIn) {
         if (nBlocks==0 || nHeightFirst > nHeightIn)
             nHeightFirst = nHeightIn;
         if (nBlocks==0 || nTimeFirst > nTimeIn)
             nTimeFirst = nTimeIn;
         nBlocks++;
         if (nHeightIn > nHeightFirst)
             nHeightLast = nHeightIn;
         if (nTimeIn > nTimeLast)
             nTimeLast = nTimeIn;
     }
};

extern CCriticalSection cs_LastBlockFile;
extern CBlockFileInfo infoLastBlockFile;
extern int nLastBlockFile;

struct CCoinsStats
{
    int nHeight;
    uint256 hashBlock;
    uint64_t nTransactions;
    uint64_t nTransactionOutputs;
    uint64_t nSerializedSize;
    uint256 hashSerialized;
    int64_t nTotalAmount;

    CCoinsStats() : nHeight(0), hashBlock(0), nTransactions(0), nTransactionOutputs(0), nSerializedSize(0), hashSerialized(0), nTotalAmount(0) {}
};



/** Abstract view on the open txout dataset. */
class CCoinsView
{
public:
    // Retrieve the CCoins (unspent transaction outputs) for a given txid
    virtual bool GetCoins(const uint256 &txid, CCoins &coins);

    // Modify the CCoins for a given txid
    virtual bool SetCoins(const uint256 &txid, const CCoins &coins);

    // Just check whether we have data for a given txid.
    // This may (but cannot always) return true for fully spent transactions
    virtual bool HaveCoins(const uint256 &txid);

    // Retrieve the block index whose state this CCoinsView currently represents
    virtual CBlockIndex *GetBestBlock();

    // Modify the currently active block index
    virtual bool SetBestBlock(CBlockIndex *pindex);

    // Do a bulk modification (multiple SetCoins + one SetBestBlock)
    virtual bool BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex);

    // Calculate statistics about the unspent transaction output set
    virtual bool GetStats(CCoinsStats &stats);

    // As we use CCoinsViews polymorphically, have a virtual destructor
    virtual ~CCoinsView() {}
};

/** CCoinsView backed by another CCoinsView */
class CCoinsViewBacked : public CCoinsView
{
protected:
    CCoinsView *base;

public:
    CCoinsViewBacked(CCoinsView &viewIn);
    bool GetCoins(const uint256 &txid, CCoins &coins);
    bool SetCoins(const uint256 &txid, const CCoins &coins);
    bool HaveCoins(const uint256 &txid);
    CBlockIndex *GetBestBlock();
    bool SetBestBlock(CBlockIndex *pindex);
    void SetBackend(CCoinsView &viewIn);
    bool BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex);
    bool GetStats(CCoinsStats &stats);
};

/** CCoinsView that adds a memory cache for transactions to another CCoinsView */
class CCoinsViewCache : public CCoinsViewBacked
{
protected:
    CBlockIndex *pindexTip;
    std::map<uint256,CCoins> cacheCoins;

public:
    CCoinsViewCache(CCoinsView &baseIn, bool fDummy = false);

    // Standard CCoinsView methods
    bool GetCoins(const uint256 &txid, CCoins &coins);
    bool SetCoins(const uint256 &txid, const CCoins &coins);
    bool HaveCoins(const uint256 &txid);
    CBlockIndex *GetBestBlock();
    bool SetBestBlock(CBlockIndex *pindex);
    bool BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex);

    // Return a modifiable reference to a CCoins. Check HaveCoins first.
    // Many methods explicitly require a CCoinsViewCache because of this method, to reduce
    // copying.
    CCoins &GetCoins(const uint256 &txid);

    // Push the modifications applied to this cache to its base.
    // Failure to call this method before destruction will cause the changes to be forgotten.
    bool Flush();

    // Calculate the size of the cache (in number of transactions)
    unsigned int GetCacheSize();

private:
    std::map<uint256,CCoins>::iterator FetchCoins(const uint256 &txid);
};

/** CCoinsView that brings transactions from a memorypool into view.
    It does not check for spendings by memory pool transactions. */
class CCoinsViewMemPool : public CCoinsViewBacked
{
protected:
    CTxMemPool &mempool;

public:
    CCoinsViewMemPool(CCoinsView &baseIn, CTxMemPool &mempoolIn);
    bool GetCoins(const uint256 &txid, CCoins &coins);
    bool HaveCoins(const uint256 &txid);
};

/** Global variable that points to the active CCoinsView (protected by cs_main) */
extern CCoinsViewCache *pcoinsTip;

/** Global variable that points to the active block tree (protected by cs_main) */
extern CBlockTreeDB *pblocktree;

struct CBlockTemplate
{
    CBlock block;
    std::vector<int64_t> vTxFees;
    std::vector<int64_t> vTxSigOps;
};

#if defined(_M_IX86) || defined(__i386__) || defined(__i386) || defined(_M_X64) || defined(__x86_64__) || defined(_M_AMD64)
extern unsigned int cpuid_edx;
#endif


/** Used to relay blocks as header + vector<merkle branch>
 * to filtered nodes.
 */
class CMerkleBlock
{
public:
    // Public only for unit testing
    CBlockHeader header;
    CPartialMerkleTree txn;

public:
    // Public only for unit testing and relay testing
    // (not relayed)
    std::vector<std::pair<unsigned int, uint256> > vMatchedTxn;

    // Create from a CBlock, filtering transactions according to filter
    // Note that this will call IsRelevantAndUpdate on the filter for each transaction,
    // thus the filter will likely be modified.
    CMerkleBlock(const CBlock& block);

    IMPLEMENT_SERIALIZE
    (
        READWRITE(header);
        READWRITE(txn);
    )
};


class CWalletInterface {
protected:
    virtual void SyncTransaction(const CTransaction &tx, const CBlock *pblock, bool fConnect) =0;
    virtual void StakeTransaction(const CScript& script, int64_t nStakeReward, bool fConnect) =0;
    virtual void EraseFromWallet(const uint256 &hash) =0;
    virtual void SetBestChain(const CBlockLocator &locator) =0;
    virtual void UpdatedTransaction(const uint256 &hash) =0;
    virtual void Inventory(const uint256 &hash) =0;
    virtual void ResendWalletTransactions(bool fForce) =0;
    friend void ::RegisterWallet(CWalletInterface*);
    friend void ::UnregisterWallet(CWalletInterface*);
    friend void ::UnregisterAllWallets();
};

#endif
