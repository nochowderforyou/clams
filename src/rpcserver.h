// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef _BITCOINRPC_SERVER_H_
#define _BITCOINRPC_SERVER_H_ 1

#include "uint256.h"
#include "rpcprotocol.h"

#include <list>
#include <map>

#include <boost/function.hpp>

#include <univalue.h>

class CBlockIndex;

void StartRPCThreads();
void StopRPCThreads();

/*
  Type-check arguments; throws JSONRPCError if wrong type given. Does not check that
  the right number of arguments are passed, just that any passed are the correct type.
  Use like:  RPCTypeCheck(params, boost::assign::list_of(str_type)(int_type)(obj_type));
*/
void RPCTypeCheck(const UniValue& params,
                  const std::list<UniValue::VType>& typesExpected, bool fAllowNull=false);
/*
  Check for expected keys/value types in an Object.
  Use like: RPCTypeCheck(object, boost::assign::map_list_of("name", str_type)("value", int_type));
*/
void RPCTypeCheckObj(const UniValue& o,
                  const std::map<std::string, UniValue::VType>& typesExpected, bool fAllowNull=false);

/*
  Run func nSeconds from now. Uses boost deadline timers.
  Overrides previous timer <name> (if any).
 */
void RPCRunLater(const std::string& name, boost::function<void(void)> func, int64_t nSeconds);


typedef UniValue(*rpcfn_type)(const UniValue& params, bool fHelp);

class CRPCCommand
{
public:
    std::string name;
    rpcfn_type actor;
    bool okSafeMode;
    bool threadSafe;
    bool reqWallet;
};

/**
 * Bitcoin RPC command dispatcher.
 */
class CRPCTable
{
private:
    std::map<std::string, const CRPCCommand*> mapCommands;
public:
    CRPCTable();
    const CRPCCommand* operator[](std::string name) const;
    std::string help(std::string name) const;

    /**
     * Execute a method.
     * @param method   Method to execute
     * @param params   Array of arguments (JSON objects)
     * @returns Result of the call.
     * @throws an exception (json_spirit::Value) when an error happens.
     */
    UniValue execute(const std::string &method, const UniValue &params) const;
};

extern const CRPCTable tableRPC;

extern void InitRPCMining();
extern void ShutdownRPCMining();

extern int64_t nWalletUnlockTime;
extern int64_t AmountFromValue(const UniValue& value);
extern UniValue ValueFromAmount(int64_t amount);
extern double GetDifficulty(const CBlockIndex* blockindex = NULL);

extern double GetPoWMHashPS();
extern double GetPoSKernelPS();

extern std::string HelpRequiringPassphrase();
extern void EnsureWalletIsUnlocked();

//
// Utilities: convert hex-encoded Values
// (throws error if not hex).
//
extern uint256 ParseHashV(const UniValue& v, std::string strName);
extern uint256 ParseHashO(const UniValue& o, std::string strKey);
extern std::vector<unsigned char> ParseHexV(const UniValue& v, std::string strName);
extern std::vector<unsigned char> ParseHexO(const UniValue& o, std::string strKey);

extern UniValue getconnectioncount(const UniValue& params, bool fHelp); // in rpcnet.cpp
extern UniValue getpeerinfo(const UniValue& params, bool fHelp);
extern UniValue ping(const UniValue& params, bool fHelp);
extern UniValue addnode(const UniValue& params, bool fHelp);
extern UniValue getaddednodeinfo(const UniValue& params, bool fHelp);
extern UniValue getnettotals(const UniValue& params, bool fHelp);

extern UniValue importwalletdump(const UniValue& params, bool fHelp); // in rpcdump.cpp
extern UniValue importprivkey(const UniValue& params, bool fHelp);
extern UniValue dumpprivkey(const UniValue& params, bool fHelp); 
extern UniValue dumpwallet(const UniValue& params, bool fHelp); 
extern UniValue importwallet(const UniValue& params, bool fHelp);
extern UniValue deleteprivkey(const UniValue& params, bool fHelp);

extern UniValue sendalert(const UniValue& params, bool fHelp);

extern UniValue getsubsidy(const UniValue& params, bool fHelp);
extern UniValue getmininginfo(const UniValue& params, bool fHelp);
extern UniValue getstakinginfo(const UniValue& params, bool fHelp);
extern UniValue getwork(const UniValue& params, bool fHelp);
extern UniValue getworkex(const UniValue& params, bool fHelp);
extern UniValue getblocktemplate(const UniValue& params, bool fHelp);
extern UniValue submitblock(const UniValue& params, bool fHelp);

extern UniValue getnewaddress(const UniValue& params, bool fHelp); // in rpcwallet.cpp
extern UniValue getaccountaddress(const UniValue& params, bool fHelp);
extern UniValue setaccount(const UniValue& params, bool fHelp);
extern UniValue getaccount(const UniValue& params, bool fHelp);
extern UniValue getaddressesbyaccount(const UniValue& params, bool fHelp);
extern UniValue sendtoaddress(const UniValue& params, bool fHelp);
extern UniValue signmessage(const UniValue& params, bool fHelp);
extern UniValue verifymessage(const UniValue& params, bool fHelp);
extern UniValue setspeech(const UniValue& params, bool fHelp);
extern UniValue setstakespeech(const UniValue& params, bool fHelp);
extern UniValue setweightedstakespeech(const UniValue& params, bool fHelp);
extern UniValue getstakedbyaddress(const UniValue& params, bool fHelp);
extern UniValue getreceivedbyaddress(const UniValue& params, bool fHelp);
extern UniValue getreceivedbyaccount(const UniValue& params, bool fHelp);
extern UniValue getbalance(const UniValue& params, bool fHelp);
extern UniValue movecmd(const UniValue& params, bool fHelp);
extern UniValue sendfrom(const UniValue& params, bool fHelp);
extern UniValue sendmany(const UniValue& params, bool fHelp);
extern UniValue addmultisigaddress(const UniValue& params, bool fHelp);
extern UniValue addredeemscript(const UniValue& params, bool fHelp);
extern UniValue listreceivedbyaddress(const UniValue& params, bool fHelp);
extern UniValue listreceivedbyaccount(const UniValue& params, bool fHelp);
extern UniValue listtransactions(const UniValue& params, bool fHelp);
extern UniValue listaddressgroupings(const UniValue& params, bool fHelp);
extern UniValue listbalances(const UniValue& params, bool fHelp);
extern UniValue listaccounts(const UniValue& params, bool fHelp);
extern UniValue listsinceblock(const UniValue& params, bool fHelp);
extern UniValue gettransaction(const UniValue& params, bool fHelp);
extern UniValue dumpbootstrap(const UniValue& params, bool fHelp);
extern UniValue backupwallet(const UniValue& params, bool fHelp);
extern UniValue keypoolrefill(const UniValue& params, bool fHelp);
extern UniValue walletpassphrase(const UniValue& params, bool fHelp);
extern UniValue walletpassphrasechange(const UniValue& params, bool fHelp);
extern UniValue walletlock(const UniValue& params, bool fHelp);
extern UniValue encryptwallet(const UniValue& params, bool fHelp);
extern UniValue validateaddress(const UniValue& params, bool fHelp);
extern UniValue validateoutputs(const UniValue& params, bool fHelp);
extern UniValue gettxoutsetinfo(const UniValue& params, bool fHelp);
extern UniValue getinfo(const UniValue& params, bool fHelp);
extern UniValue reservebalance(const UniValue& params, bool fHelp);
extern UniValue resendtx(const UniValue& params, bool fHelp);
extern UniValue makekeypair(const UniValue& params, bool fHelp);
extern UniValue validatepubkey(const UniValue& params, bool fHelp);
extern UniValue getnewpubkey(const UniValue& params, bool fHelp);

extern UniValue getrawtransaction(const UniValue& params, bool fHelp); // in rcprawtransaction.cpp
extern UniValue listunspent(const UniValue& params, bool fHelp);
extern UniValue createrawtransaction(const UniValue& params, bool fHelp);
extern UniValue decoderawtransaction(const UniValue& params, bool fHelp);
extern UniValue decodescript(const UniValue& params, bool fHelp);
extern UniValue signrawtransaction(const UniValue& params, bool fHelp);
extern UniValue sendrawtransaction(const UniValue& params, bool fHelp);

extern UniValue getbestblockhash(const UniValue& params, bool fHelp); // in rpcblockchain.cpp
extern UniValue getblockcount(const UniValue& params, bool fHelp); // in rpcblockchain.cpp
extern UniValue getdifficulty(const UniValue& params, bool fHelp);
extern UniValue settxfee(const UniValue& params, bool fHelp);
extern UniValue getrawmempool(const UniValue& params, bool fHelp);
extern UniValue getblockhash(const UniValue& params, bool fHelp);
extern UniValue getblock(const UniValue& params, bool fHelp);
extern UniValue getblockbynumber(const UniValue& params, bool fHelp);
extern UniValue getcheckpoint(const UniValue& params, bool fHelp);
extern UniValue getstaketo(const UniValue& params, bool fHelp);
extern UniValue setstaketo(const UniValue& params, bool fHelp);
extern UniValue getrewardto(const UniValue& params, bool fHelp);
extern UniValue setrewardto(const UniValue& params, bool fHelp);

extern UniValue sendnotarytransaction(const UniValue& params, bool fHelp);
extern UniValue createclamour(const UniValue& params, bool fHelp);
extern UniValue getnotarytransaction(const UniValue& params, bool fHelp);

extern UniValue getclamour(const UniValue& params, bool fHelp);
extern UniValue listclamours(const UniValue& params, bool fHelp);
extern UniValue getsupport(const UniValue& params, bool fHelp);

#endif
