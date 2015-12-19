// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2013 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "base58.h"
#include "clamspeech.h"
#include "init.h"
#include "main.h"
#include "net.h"
#include "netbase.h"
#include "rpcserver.h"
#include "txdb.h"
#include "timedata.h"
#include "util.h"
#ifdef ENABLE_WALLET
#include "wallet.h"
#include "walletdb.h"
#endif

#include <stdint.h>
#include <fstream>

#include <boost/assign/list_of.hpp>
#include <univalue.h>

using namespace std;
using namespace boost;
using namespace boost::assign;

typedef map<string, CClamour*> mapClamour_t;
extern  mapClamour_t mapClamour;

UniValue getinfo(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getinfo\n"
            "Returns an object containing various state info.");

    proxyType proxy;
    GetProxy(NET_IPV4, proxy);

    UniValue obj(UniValue::VOBJ);
    UniValue diff(UniValue::VOBJ);
    obj.push_back(Pair("version",       FormatFullVersion()));
    obj.push_back(Pair("protocolversion",(int)PROTOCOL_VERSION));
#ifdef ENABLE_WALLET
    if (pwalletMain) {
        obj.push_back(Pair("walletversion", pwalletMain->GetVersion()));
        obj.push_back(Pair("balance",       ValueFromAmount(pwalletMain->GetBalance())));
        obj.push_back(Pair("newmint",       ValueFromAmount(pwalletMain->GetNewMint())));
        obj.push_back(Pair("stake",         ValueFromAmount(pwalletMain->GetStake())));
    }
#endif
    obj.push_back(Pair("blocks",        (int)nBestHeight));
    obj.push_back(Pair("timeoffset",    (int64_t)GetTimeOffset()));
    obj.push_back(Pair("moneysupply",   ValueFromAmount(pindexBest->nMoneySupply)));
    obj.push_back(Pair("digsupply",     ValueFromAmount(pindexBest->nDigsupply)));
    obj.push_back(Pair("stakesupply",   ValueFromAmount(pindexBest->nStakeSupply)));
    obj.push_back(Pair("activesupply",  ValueFromAmount(pindexBest->nDigsupply + pindexBest->nStakeSupply)));
    obj.push_back(Pair("connections",   (int)vNodes.size()));
    obj.push_back(Pair("proxy",         (proxy.first.IsValid() ? proxy.first.ToStringIPPort() : string())));
    obj.push_back(Pair("ip",            GetLocalAddress(NULL).ToStringIP()));

    diff.push_back(Pair("proof-of-stake", GetDifficulty(GetLastBlockIndex(pindexBest, true))));
    obj.push_back(Pair("difficulty",    diff));

    obj.push_back(Pair("testnet",       TestNet()));
#ifdef ENABLE_WALLET
    if (pwalletMain) {
        obj.push_back(Pair("keypoololdest", (int64_t)pwalletMain->GetOldestKeyPoolTime()));
        obj.push_back(Pair("keypoolsize",   (int)pwalletMain->GetKeyPoolSize()));
    }
    obj.push_back(Pair("paytxfee",      ValueFromAmount(nTransactionFee)));
    obj.push_back(Pair("mininput",      ValueFromAmount(nMinimumInputValue)));
    if (pwalletMain && pwalletMain->IsCrypted())
        obj.push_back(Pair("unlocked_until", (int64_t)nWalletUnlockTime));
#endif
    obj.push_back(Pair("errors",        GetWarnings("statusbar")));
    return obj;
}

#ifdef ENABLE_WALLET
class DescribeAddressVisitor : public boost::static_visitor<UniValue>
{
public:
    UniValue operator()(const CNoDestination &dest) const { return UniValue(); }

    UniValue operator()(const CKeyID &keyID) const {
        UniValue obj(UniValue::VOBJ);
        CPubKey vchPubKey;
        pwalletMain->GetPubKey(keyID, vchPubKey);
        obj.push_back(Pair("isscript", false));
        obj.push_back(Pair("pubkey", HexStr(vchPubKey)));
        obj.push_back(Pair("iscompressed", vchPubKey.IsCompressed()));
        return obj;
    }

    UniValue operator()(const CScriptID &scriptID) const {
        UniValue obj(UniValue::VOBJ);
        obj.push_back(Pair("isscript", true));
        CScript subscript;
        pwalletMain->GetCScript(scriptID, subscript);
        std::vector<CTxDestination> addresses;
        txnouttype whichType;
        int nRequired;
        ExtractDestinations(subscript, whichType, addresses, nRequired);
        obj.push_back(Pair("script", GetTxnOutputType(whichType)));
        obj.push_back(Pair("hex", HexStr(subscript.begin(), subscript.end())));
        UniValue a(UniValue::VARR);
        BOOST_FOREACH(const CTxDestination& addr, addresses)
            a.push_back(CBitcoinAddress(addr).ToString());
        obj.push_back(Pair("addresses", a));
        if (whichType == TX_MULTISIG)
            obj.push_back(Pair("sigsrequired", nRequired));
        return obj;
    }
};
#endif

UniValue validateaddress(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "validateaddress <clamaddress>\n"
            "Return information about <clamaddress>.");

    CBitcoinAddress address(params[0].get_str());
    bool isValid = address.IsValid();

    UniValue ret(UniValue::VOBJ);
    ret.push_back(Pair("isvalid", isValid));
    if (isValid)
    {
        CTxDestination dest = address.Get();
        string currentAddress = address.ToString();
        ret.push_back(Pair("address", currentAddress));
#ifdef ENABLE_WALLET
        bool fMine = pwalletMain ? IsMine(*pwalletMain, dest) : false;
        ret.push_back(Pair("ismine", fMine));
        if (fMine) {
            UniValue detail = boost::apply_visitor(DescribeAddressVisitor(), dest);
            ret.pushKVs(detail);
        }
        if (pwalletMain && pwalletMain->mapAddressBook.count(dest))
            ret.push_back(Pair("account", pwalletMain->mapAddressBook[dest]));
#endif
    }
    return ret;
}


static void validateoutputs_check_unconfirmed_spend(COutPoint& outpoint, CTransaction& tx, UniValue& entry)
{
    // check whether unconfirmed output is already spent
    LOCK(mempool.cs); // protect mempool.mapNextTx
    if (mempool.mapNextTx.count(outpoint)) {
        // pull details from mempool
        CInPoint in = mempool.mapNextTx[outpoint];
        UniValue details(UniValue::VOBJ);
        entry.push_back(Pair("status", "spent"));
        details.push_back(Pair("txid", in.ptx->GetHash().GetHex()));
        details.push_back(Pair("vin", int(in.n)));
        details.push_back(Pair("confirmations", 0));
        entry.push_back(Pair("spent", details));
    } else {
        entry.push_back(Pair("status", "unspent"));
    }
}

UniValue validateoutputs(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "validateoutputs [{\"txid\":txid,\"vout\":n},...]\n"
            "Return information about outputs (whether they exist, and whether they have been spent).");

    UniValue inputs = params[0].get_array();

    UniValue ret(UniValue::VARR);

    for (unsigned int idx = 0; idx < inputs.size(); idx++)
    //BOOST_FOREACH(UniValue::VOBJ input, inputs)
    {
        UniValue o = inputs[idx].get_obj();

        const UniValue& txid_v = find_value(o, "txid");
        if (txid_v.isNull())
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, missing txid key");
        string txid = txid_v.get_str();
        if (!IsHex(txid))
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, expected hex txid");

        const UniValue& vout_v = find_value(o, "vout");
        if (!vout_v.isNum())
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, missing vout key");
        int nOutput = vout_v.get_int();
        if (nOutput < 0)
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, vout must be positive");

        UniValue entry(UniValue::VOBJ);
        entry.push_back(Pair("txid", txid));
        entry.push_back(Pair("vout", nOutput));

        CTransaction tx;
        CCoins coins;
        uint256 hashBlock = 0;
        COutPoint outpoint(uint256(txid), nOutput);

        if (!pcoinsTip->GetCoins(uint256(txid), coins)) {
             entry.push_back(Pair("status", "spent"));
             ret.push_back(entry);
             continue;
        }

        if (!GetTransaction(uint256(txid), tx, hashBlock, true)) {
            entry.push_back(Pair("status", "txid not found on disk"));
            ret.push_back(entry);
            continue;
        }

        // check that the output number is in range
        if ((unsigned)nOutput >= coins.vout.size()) {
            entry.push_back(Pair("status", "vout too high"));
            entry.push_back(Pair("outputs", (int)coins.vout.size()));
            ret.push_back(entry);
            continue;
        }

        entry.push_back(Pair("amount", ValueFromAmount(tx.vout[nOutput].nValue)));

        // get the address and account
        CTxDestination address;
        if (ExtractDestination(tx.vout[nOutput].scriptPubKey, address))
        {
            entry.push_back(Pair("address", CBitcoinAddress(address).ToString()));
            if (pwalletMain->mapAddressBook.count(address))
                entry.push_back(Pair("account", pwalletMain->mapAddressBook[address]));
        }

        const CScript& pk = tx.vout[outpoint.n].scriptPubKey;
        entry.push_back(Pair("scriptPubKey", HexStr(pk.begin(), pk.end())));

        // is the output confirmed?
        if (hashBlock == 0) {
            entry.push_back(Pair("confirmations", 0));
            validateoutputs_check_unconfirmed_spend(outpoint, tx, entry);
            ret.push_back(entry);
            continue;
        }

        // find the block containing the output
        map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashBlock);
        if (mi != mapBlockIndex.end() && (*mi).second)
        {
            CBlockIndex* pindex = (*mi).second;
            if (pindex->IsInMainChain()) {
                entry.push_back(Pair("height", pindex->nHeight));
                entry.push_back(Pair("confirmations", 1 + nBestHeight - pindex->nHeight));
            } else {
                LogPrintf("can't find block with hash %s\n", hashBlock.GetHex().c_str());
                entry.push_back(Pair("confirmations", 0));
            }
        }

        // check whether any confirmed transaction spends this output
        if (!coins.vout[nOutput].IsNull()) {
            // if not, check for an unconfirmed spend
            validateoutputs_check_unconfirmed_spend(outpoint, tx, entry);
            ret.push_back(entry);
            continue;
        }

        entry.push_back(Pair("status", "spent"));
        
        ret.push_back(entry);
    }
    
    return ret;
}


UniValue validatepubkey(const UniValue& params, bool fHelp)
{
    if (fHelp || !params.size() || params.size() > 2)
        throw runtime_error(
            "validatepubkey <clampubkey>\n"
            "Return information about <clampubkey>.");

    std::vector<unsigned char> vchPubKey = ParseHex(params[0].get_str());
    CPubKey pubKey(vchPubKey);

    bool isValid = pubKey.IsValid();
    bool isCompressed = pubKey.IsCompressed();
    CKeyID keyID = pubKey.GetID();

    CBitcoinAddress address;
    address.Set(keyID);

    UniValue ret(UniValue::VOBJ);
    ret.push_back(Pair("isvalid", isValid));
    if (isValid)
    {
        CTxDestination dest = address.Get();
        string currentAddress = address.ToString();
        ret.push_back(Pair("address", currentAddress));
        ret.push_back(Pair("iscompressed", isCompressed));
#ifdef ENABLE_WALLET
        bool fMine = pwalletMain ? IsMine(*pwalletMain, dest) : false;
        ret.push_back(Pair("ismine", fMine));
        if (fMine) {
            UniValue detail = boost::apply_visitor(DescribeAddressVisitor(), dest);
            ret.pushKVs(detail);
        }
        if (pwalletMain && pwalletMain->mapAddressBook.count(dest))
            ret.push_back(Pair("account", pwalletMain->mapAddressBook[dest]));
#endif
    }
    return ret;
}

UniValue verifymessage(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 3)
        throw runtime_error(
            "verifymessage <clamaddress> <signature> <message>\n"
            "Verify a signed message");

    string strAddress  = params[0].get_str();
    string strSign     = params[1].get_str();
    string strMessage  = params[2].get_str();

    CBitcoinAddress addr(strAddress);
    if (!addr.IsValid())
        throw JSONRPCError(RPC_TYPE_ERROR, "Invalid address");

    CKeyID keyID;
    if (!addr.GetKeyID(keyID))
        throw JSONRPCError(RPC_TYPE_ERROR, "Address does not refer to key");

    bool fInvalid = false;
    vector<unsigned char> vchSig = DecodeBase64(strSign.c_str(), &fInvalid);

    if (fInvalid)
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Malformed base64 encoding");

    CHashWriter ss(SER_GETHASH, 0);
    ss << strMessageMagic;
    ss << strMessage;

    CPubKey pubkey;
    if (!pubkey.RecoverCompact(ss.GetHash(), vchSig))
        return false;

    return (pubkey.GetID() == keyID);
}

UniValue setspeech(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "setspeech <text>\n"
            "Sets the text to be used as the transaction comment when making transactions.");

    strDefaultSpeech = params[0].get_str();

    LogPrint("speech", "set default speech to \"%s\"\n", strDefaultSpeech);

    return NullUniValue;
}

UniValue setstakespeech(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "setstakespeech <text>\n"
            "Sets the text to be as the transaction comment when staking");

    strDefaultStakeSpeech = params[0].get_str();

    LogPrint("stakespeech", "set default stakespeech to \"%s\"\n", strDefaultStakeSpeech);

    return NullUniValue;
}

UniValue setweightedstakespeech(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() > 1)
        throw runtime_error(
            "setweightedstakespeech [path]\n"
            "Loads a file containing a list of texts to be as the transaction comment when staking.\n"
            "Each line in the file should contain a positive integer specifying the probabalistic weight for that line, then a space, then the speech.\n"
            "If no path is provided or any errors occur opening or parsing the file then weighted staking isn't used at all.");

    weightedStakeSpeech.clear();

    if (params.size() == 0)
        return NullUniValue;

    string strPath = params[0].get_str();

    if (!filesystem::exists(strPath))
        throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter: file doesn't exist");

    std::ifstream speechfile(strPath.c_str());

    if(!speechfile) //Always test the file open.
        throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter: can't open file");

    string line;
    size_t pos;
    const char *start;
    char *end;
    int count = 0;
    while (getline(speechfile, line, '\n')) {
        count++;
        start = line.c_str();

        // blank line and lines starting with '#' are comments
        if (*start == '\0' || *start == '#')
            continue;

        unsigned long weight = strtoul(start, &end, 10);
        if (weight == ULONG_MAX && errno == ERANGE) {
            weightedStakeSpeech.clear();
            throw JSONRPCError(RPC_INVALID_PARAMETER, strprintf("Weight out of range, line %d", count));
        }

        pos = end - start;
        if (pos == 0) {
            weightedStakeSpeech.clear();
            throw JSONRPCError(RPC_INVALID_PARAMETER, strprintf("Invalid weight, line %d", count));
        }

        if (*end == '\0')
            line = string();
        else if (*end == ' ')
            line = line.substr(pos+1);
        else {
            weightedStakeSpeech.clear();
            throw JSONRPCError(RPC_INVALID_PARAMETER, strprintf("No space after weight, line %d", count));
        }

        weightedStakeSpeech.insert(weight, line);
    }

    return strprintf("loaded %d weighted stake speech text(s)", weightedStakeSpeech.size());
}

UniValue getclamour(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "getclamour <pid>\n"
            "Returns an object containing info about the specified petition ID");

    string pid = params[0].get_str();

    map<string, CClamour*>::iterator mi = mapClamour.find(pid);
    if (mi == mapClamour.end())
        return NullUniValue;

    UniValue ret(UniValue::VOBJ);
    CClamour *clamour = mi->second;

    ret.push_back(Pair("pid", pid));
    ret.push_back(Pair("hash", clamour->strHash));
    if (clamour->strURL.length())
        ret.push_back(Pair("url", clamour->strURL));
    ret.push_back(Pair("txid", clamour->txid.GetHex()));
    ret.push_back(Pair("confirmations", pindexBest->nHeight - clamour->nHeight + 1));
    return ret;
}

UniValue listclamours(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() > 2)
        throw runtime_error(
            "listclamours [minconf=1] [maxconf=9999999]\n"
            "Returns an array of objects containing info about all registered petitions\n"
            "with between minconf and maxconf (inclusive) confirmations.");

    RPCTypeCheck(params, list_of(UniValue::VNUM)(UniValue::VNUM));

    int nMinDepth = 1;
    if (params.size() > 0)
        nMinDepth = params[0].get_int();

    int nMaxDepth = 9999999;
    if (params.size() > 1)
        nMaxDepth = params[1].get_int();

    UniValue ret(UniValue::VARR);

    BOOST_FOREACH(const mapClamour_t::value_type pair, mapClamour)
    {
        CClamour *clamour = pair.second;
        int nDepth = pindexBest->nHeight - clamour->nHeight + 1;

        if (nDepth < nMinDepth || nDepth > nMaxDepth)
            continue;

        UniValue entry(UniValue::VOBJ);

        entry.push_back(Pair("pid", clamour->strHash.substr(0, 8)));
        entry.push_back(Pair("hash", clamour->strHash));
        if (clamour->strURL.length())
            entry.push_back(Pair("url", clamour->strURL));
        entry.push_back(Pair("txid", clamour->txid.GetHex()));
        entry.push_back(Pair("confirmations", nDepth));

        ret.push_back(entry);
    }

    return ret;
}

UniValue getsupport(const UniValue& params, bool fHelp)
{
    if (fHelp || params.size() > 3)
        throw runtime_error(
            "getsupport [threshold=0] [window=10000] [block=<bestblock>]\n"
            "Returns an object detailing the number of blocks supporting CLAMour petitions\n"
            "<threshold> sets a percentage threshold of support below which petitions are ignored.\n"
            "<window> sets the number of blocks to count and defaults to 10000.\n"
            "<block> sets which block ends the window, and defaults to the last block in the chain.");

    RPCTypeCheck(params, list_of(UniValue::VNUM)(UniValue::VNUM)(UniValue::VNUM));

    double dThreshold;
    int nWindow, nBlock;
    map<string,int> mapSupport;
    typedef pair<string,int> mapSupport_pair;

    if (params.size() > 0) {
        dThreshold = params[0].get_real();
        if (dThreshold < 0 || dThreshold > 100)
            throw runtime_error("Threshold percentage out of range.");
    } else
        dThreshold = 0;

    if (params.size() > 1)
        nWindow = params[1].get_int();
    else
        nWindow = 10000;

    if (params.size() > 2) {
        nBlock = params[2].get_int();
        if (nBlock < 0 || nBlock > nBestHeight)
            throw runtime_error("Block number out of range.");
    } else
        nBlock = nBestHeight;

    if (nWindow < 1)
        throw runtime_error("Window size must be at least 1.");

    if (nWindow > nBlock + 1)
        throw runtime_error("Window starts before block 0.");

    for (int i = nBlock + 1 - nWindow; i <= nBlock; i++) {
        CBlockIndex* pblockindex = FindBlockByHeight(i);
        set<string> sup = pblockindex->GetSupport();
        BOOST_FOREACH(const string &s, sup) {
            // LogPrintf("%d supports %s\n", i, s);
            mapSupport[s]++;
        }
    }

    UniValue ret(UniValue::VOBJ);
    UniValue counts(UniValue::VOBJ);
    ret.push_back(Pair("threshold", dThreshold));
    ret.push_back(Pair("window", nWindow));
    ret.push_back(Pair("endblock", nBlock));
    ret.push_back(Pair("startblock", nBlock + 1 - nWindow));
    BOOST_FOREACH(const mapSupport_pair &p, mapSupport)
        if (p.second * 100 >= dThreshold * nWindow)
            counts.push_back(Pair(p.first, p.second));
    ret.push_back(Pair("support", counts));

    return ret;
}
