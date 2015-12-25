// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "rpcserver.h"
#include "main.h"
#include "kernel.h"
#include "checkpoints.h"

#include <univalue.h>

using namespace std;

extern void TxToJSON(const CTransaction& tx, const uint256 hashBlock, UniValue& entry);
extern enum Checkpoints::CPMode CheckpointsMode;

double GetDifficulty(const CBlockIndex* blockindex)
{
    // Floating point number that is a multiple of the minimum difficulty,
    // minimum difficulty = 1.0.
    if (blockindex == NULL)
    {
        if (pindexBest == NULL)
            return 1.0;
        else
            blockindex = GetLastBlockIndex(pindexBest, false);
    }

    int nShift = (blockindex->nBits >> 24) & 0xff;

    double dDiff =
        (double)0x0000ffff / (double)(blockindex->nBits & 0x00ffffff);

    while (nShift < 29)
    {
        dDiff *= 256.0;
        nShift++;
    }
    while (nShift > 29)
    {
        dDiff /= 256.0;
        nShift--;
    }

    return dDiff;
}

double GetPoWMHashPS()
{
    if (pindexBest->nHeight >= LAST_POW_BLOCK)
        return 0;

    int nPoWInterval = 72;
    int64_t nTargetSpacingWorkMin = 30, nTargetSpacingWork = 30;

    CBlockIndex* pindex = pindexGenesisBlock;
    CBlockIndex* pindexPrevWork = pindexGenesisBlock;

    while (pindex)
    {
        if (pindex->IsProofOfWork())
        {
            int64_t nActualSpacingWork = pindex->GetBlockTime() - pindexPrevWork->GetBlockTime();
            nTargetSpacingWork = ((nPoWInterval - 1) * nTargetSpacingWork + nActualSpacingWork + nActualSpacingWork) / (nPoWInterval + 1);
            nTargetSpacingWork = max(nTargetSpacingWork, nTargetSpacingWorkMin);
            pindexPrevWork = pindex;
        }

        pindex = pindex->pnext;
    }

    return GetDifficulty() * 4294.967296 / nTargetSpacingWork;
}

double GetPoSKernelPS()
{
    int nPoSInterval = 72;
    double dStakeKernelsTriedAvg = 0;
    int nStakesHandled = 0, nStakesTime = 0;

    CBlockIndex* pindex = pindexBest;
    CBlockIndex* pindexPrevStake = NULL;

    while (pindex && nStakesHandled <= nPoSInterval)
    {
        if (pindex->IsProofOfStake())
        {
            if (nStakesHandled < nPoSInterval)
                dStakeKernelsTriedAvg += GetDifficulty(pindex) * 0x100000000;
            nStakesTime += pindexPrevStake ? (pindexPrevStake->nTime - pindex->nTime) : 0;
            pindexPrevStake = pindex;
            nStakesHandled++;
        }

        pindex = pindex->pprev;
    }

    double result = 0;

    if (nStakesTime)
        result = dStakeKernelsTriedAvg / nStakesTime;

    if (IsProtocolV2(nBestHeight))
        result *= STAKE_TIMESTAMP_MASK + 1;

    return result;
}

UniValue blockToJSON(const CBlock& block, const CBlockIndex* blockindex, bool fPrintTransactionDetail)
{
    UniValue result(UniValue::VOBJ);
    result.push_back(Pair("hash", block.GetHash().GetHex()));
    int confirmations = -1;
    // Only report confirmations if the block is on the main chain
    if (blockindex->IsInMainChain())
        confirmations = nBestHeight - blockindex->nHeight + 1;
    result.push_back(Pair("confirmations", confirmations));
    result.push_back(Pair("size", (int)::GetSerializeSize(block, SER_NETWORK, PROTOCOL_VERSION)));
    result.push_back(Pair("height", blockindex->nHeight));
    result.push_back(Pair("version", block.nVersion));
    result.push_back(Pair("merkleroot", block.hashMerkleRoot.GetHex()));
    result.push_back(Pair("mint", ValueFromAmount(blockindex->nMint)));
    result.push_back(Pair("moneysupply", ValueFromAmount(blockindex->nMoneySupply)));
    result.push_back(Pair("digsupply", ValueFromAmount(blockindex->nDigsupply)));
    result.push_back(Pair("stakesupply", ValueFromAmount(blockindex->nStakeSupply)));
    result.push_back(Pair("activesupply", ValueFromAmount(blockindex->nDigsupply + blockindex->nStakeSupply)));
    result.push_back(Pair("time", (int64_t)block.GetBlockTime()));
    result.push_back(Pair("nonce", (uint64_t)block.nNonce));
    result.push_back(Pair("bits", strprintf("%08x", block.nBits)));
    result.push_back(Pair("difficulty", GetDifficulty(blockindex)));
    result.push_back(Pair("blocktrust", leftTrim(blockindex->GetBlockWork().getuint256().ToString(), '0')));
    result.push_back(Pair("chaintrust", leftTrim(blockindex->nChainWork.GetHex(), '0')));
    if (blockindex->pprev)
        result.push_back(Pair("previousblockhash", blockindex->pprev->GetBlockHash().GetHex()));
    if (blockindex->pnext)
        result.push_back(Pair("nextblockhash", blockindex->pnext->GetBlockHash().GetHex()));

    result.push_back(Pair("flags", strprintf("%s%s", blockindex->IsProofOfStake()? "proof-of-stake" : "proof-of-work", blockindex->GeneratedStakeModifier()? " stake-modifier": "")));
    result.push_back(Pair("proofhash", blockindex->hashProof.GetHex()));
    result.push_back(Pair("entropybit", (int)blockindex->GetStakeEntropyBit()));
    result.push_back(Pair("modifier", strprintf("%016x", blockindex->nStakeModifier)));
    UniValue clamours(UniValue::VARR);
    BOOST_FOREACH (const CClamour& clamour, blockindex->vClamour)
    {
        UniValue entry(UniValue::VOBJ);
        entry.push_back(Pair("txid", clamour.txid.GetHex()));
        entry.push_back(Pair("hash", clamour.strHash));
        if (clamour.strURL.length())
            entry.push_back(Pair("url", clamour.strURL));
        clamours.push_back(entry);
    }
    if (clamours.size())
        result.push_back(Pair("clamours", clamours));

    UniValue txinfo(UniValue::VARR);
    BOOST_FOREACH (const CTransaction& tx, block.vtx)
    {
        if (fPrintTransactionDetail)
        {
            UniValue entry(UniValue::VOBJ);

            entry.push_back(Pair("txid", tx.GetHash().GetHex()));
            TxToJSON(tx, 0, entry);

            txinfo.push_back(entry);
        }
        else
            txinfo.push_back(tx.GetHash().GetHex());
    }

    result.push_back(Pair("tx", txinfo));

    if (block.IsProofOfStake())
        result.push_back(Pair("signature", HexStr(block.vchBlockSig.begin(), block.vchBlockSig.end())));

    return result;
}

UniValue ignorenextblock(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "ignorenextblock"
            );

    nBlocksToIgnore++;

    return "";
}

UniValue getbestblockhash(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getbestblockhash\n"
            "Returns the hash of the best block in the longest block chain.");

    return hashBestChain.GetHex();
}

UniValue dumpbootstrap(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() < 2 || params.size() > 3)
        throw runtime_error(
            "dumpbootstrap <destination> <endblock> [startblock=0]\n"
            "Creates a bootstrap format block dump of the blockchain in destination, which can be a directory or a path with filename, up to the given endblock number.\n"
            "Optional <startblock> is the first block number to dump.");

    string strDest = params[0].get_str();
    int nEndBlock = params[1].get_int();
    if (nEndBlock < 0 || nEndBlock > nBestHeight)
        throw runtime_error("End block number out of range.");
    int nStartBlock = 0;
    if (params.size() > 2)
        nStartBlock = params[2].get_int();
    if (nStartBlock < 0 || nStartBlock > nEndBlock)
        throw runtime_error("Start block number out of range.");

    boost::filesystem::path pathDest(strDest);
    if (boost::filesystem::is_directory(pathDest))
        pathDest /= "bootstrap.dat";

    try {
        FILE* file = fopen(pathDest.string().c_str(), "wb");
        if (!file)
            throw JSONRPCError(RPC_MISC_ERROR, "Error: Could not open bootstrap file for writing.");

        CAutoFile fileout = CAutoFile(file, SER_DISK, CLIENT_VERSION);
        if (!fileout)
            throw JSONRPCError(RPC_MISC_ERROR, "Error: Could not open bootstrap file for writing.");

        for (int nHeight = nStartBlock; nHeight <= nEndBlock; nHeight++)
        {
            CBlock block;
            CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
            block.ReadFromDisk(pblockindex);
            fileout << FLATDATA(Params().MessageStart()) << fileout.GetSerializeSize(block) << block;
        }
    } catch(const boost::filesystem::filesystem_error &e) {
        throw JSONRPCError(RPC_MISC_ERROR, "Error: Bootstrap dump failed!");
    }

    return NullUniValue;
}

UniValue getblockcount(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getblockcount\n"
            "Returns the number of blocks in the longest block chain.");

    return nBestHeight;
}


UniValue getdifficulty(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getdifficulty\n"
            "Returns the difficulty as a multiple of the minimum difficulty.");

    UniValue obj(UniValue::VOBJ);
    obj.push_back(Pair("proof-of-work",        GetDifficulty()));
    obj.push_back(Pair("proof-of-stake",       GetDifficulty(GetLastBlockIndex(pindexBest, true))));
    return obj;
}


UniValue getrawmempool(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getrawmempool\n"
            "Returns all transaction ids in memory pool.");

    vector<uint256> vtxid;
    mempool.queryHashes(vtxid);

    UniValue a(UniValue::VARR);
    BOOST_FOREACH(const uint256& hash, vtxid)
        a.push_back(hash.ToString());

    return a;
}

UniValue getblockhash(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "getblockhash <index>\n"
            "Returns hash of block in best-block-chain at <index>.");

    int nHeight = params[0].get_int();
    if (nHeight < 0 || nHeight > nBestHeight)
        throw runtime_error("Block number out of range.");

    CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
    return pblockindex->phashBlock->GetHex();
}

UniValue getblock(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 3)
        throw runtime_error(
            "getblock <hash> [txinfo] [raw]\n"
            "txinfo optional to print more detailed tx info\n"
            "raw optional to return block as raw hex data\n"
            "Returns details of a block with given block-hash.");

    std::string strHash = params[0].get_str();
    size_t len = strHash.length();

    if (len < 8)
        throw runtime_error("Please provide at least the first 8 digits of the block hash");

    CBlockIndex* pblockindex;

    if (len == 64) {
        uint256 hash(strHash);

        if (mapBlockIndex.count(hash) == 0)
            throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Block not found");

        pblockindex = mapBlockIndex[hash];
    } else {
        pblockindex = NULL;

        // std::map<uint256, CBlockIndex*> mapBlockIndex;
        BOOST_FOREACH(const PAIRTYPE(uint256, CBlockIndex*)& item, mapBlockIndex) {
            if (item.first.GetHex().substr(0, len) == strHash) {
                pblockindex = item.second;
                break;
            }
        }

        if (!pblockindex)
            throw runtime_error("Block hash not found");
    }

    CBlock block;
    block.ReadFromDisk(pblockindex);

    if (params.size() > 2 && params[2].get_bool())
    {
        CDataStream ssBlock(SER_NETWORK, PROTOCOL_VERSION);
        ssBlock << block;
        return HexStr(ssBlock.begin(), ssBlock.end());
    }

    return blockToJSON(block, pblockindex, params.size() > 1 ? params[1].get_bool() : false);
}

UniValue getblockbynumber(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 3)
        throw runtime_error(
		     "getblockbynumber <number> [txinfo] [raw]\n"
		     "txinfo optional to print more detailed tx info\n"
		     "raw optional to return block as raw hex data\n"
		     "Returns details of a block with given block-number.");

    int nHeight = params[0].get_int();
    if (nHeight < 0 || nHeight > nBestHeight)
        throw runtime_error("Block number out of range.");

    CBlock block;
    CBlockIndex* pblockindex = mapBlockIndex[hashBestChain];
    while (pblockindex->nHeight > nHeight)
        pblockindex = pblockindex->pprev;

    uint256 hash = *pblockindex->phashBlock;

    pblockindex = mapBlockIndex[hash];
    block.ReadFromDisk(pblockindex);

    if (params.size() > 2 && params[2].get_bool())
    {
        CDataStream ssBlock(SER_NETWORK, PROTOCOL_VERSION);
        ssBlock << block;
        return HexStr(ssBlock.begin(), ssBlock.end());
    }

    return blockToJSON(block, pblockindex, params.size() > 1 ? params[1].get_bool() : false);
}

// ppcoin: get information of sync-checkpoint
UniValue getcheckpoint(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getcheckpoint\n"
            "Show info of synchronized checkpoint.\n");

    UniValue result(UniValue::VOBJ);
    CBlockIndex* pindexCheckpoint;

    result.push_back(Pair("synccheckpoint", Checkpoints::hashSyncCheckpoint.ToString().c_str()));
    pindexCheckpoint = mapBlockIndex[Checkpoints::hashSyncCheckpoint];
    result.push_back(Pair("height", pindexCheckpoint->nHeight));
    result.push_back(Pair("timestamp", DateTimeStrFormat(pindexCheckpoint->GetBlockTime()).c_str()));

    // Check that the block satisfies synchronized checkpoint
    if (CheckpointsMode == Checkpoints::STRICT)
        result.push_back(Pair("policy", "strict"));

    if (CheckpointsMode == Checkpoints::ADVISORY)
        result.push_back(Pair("policy", "advisory"));

    if (CheckpointsMode == Checkpoints::PERMISSIVE)
        result.push_back(Pair("policy", "permissive"));

    if (mapArgs.count("-checkpointkey"))
        result.push_back(Pair("checkpointmaster", true));

    return result;
}

UniValue gettxoutsetinfo(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "gettxoutsetinfo\n"
            "Returns statistics about the unspent transaction output set.");

    UniValue ret(UniValue::VOBJ);

    CCoinsStats stats;
    if (pcoinsTip->GetStats(stats)) {
        ret.push_back(Pair("height", (boost::int64_t)stats.nHeight));
        ret.push_back(Pair("bestblock", stats.hashBlock.GetHex()));
        ret.push_back(Pair("transactions", (boost::int64_t)stats.nTransactions));
        ret.push_back(Pair("txouts", (boost::int64_t)stats.nTransactionOutputs));
        ret.push_back(Pair("bytes_serialized", (boost::int64_t)stats.nSerializedSize));
        ret.push_back(Pair("hash_serialized", stats.hashSerialized.GetHex()));
        ret.push_back(Pair("total_amount", ValueFromAmount(stats.nTotalAmount)));
    }
    return ret;
}
