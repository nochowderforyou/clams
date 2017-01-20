// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2013 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include <set>
#include "rpcclient.h"

#include "rpcprotocol.h"
#include "util.h"
#include "ui_interface.h"
#include "chainparams.h" // for Params().RPCPort()

#include <stdint.h>

#include <boost/asio.hpp>
#include <boost/asio/ssl.hpp>
#include <boost/algorithm/string/case_conv.hpp> 
#include <univalue.h>

using namespace std;
using namespace boost;
using namespace boost::asio;

UniValue CallRPC(const string& strMethod, const UniValue& params)
{
    if (mapArgs["-rpcuser"] == "" && mapArgs["-rpcpassword"] == "")
        throw runtime_error(strprintf(
            _("You must set rpcpassword=<password> in the configuration file:\n%s\n"
              "If the file does not exist, create it with owner-readable-only file permissions."),
                GetConfigFile().string()));

    // Connect to localhost
    bool fUseSSL = GetBoolArg("-rpcssl", false);
    asio::io_service io_service;
    ssl::context context(io_service, ssl::context::sslv23);
    context.set_options(ssl::context::no_sslv2);
    asio::ssl::stream<asio::ip::tcp::socket> sslStream(io_service, context);
    SSLIOStreamDevice<asio::ip::tcp> d(sslStream, fUseSSL);
    iostreams::stream< SSLIOStreamDevice<asio::ip::tcp> > stream(d);

    bool fWait = GetBoolArg("-rpcwait", false); // -rpcwait means try until server has started
    do {
        bool fConnected = d.connect(GetArg("-rpcconnect", "127.0.0.1"), GetArg("-rpcport", itostr(Params().RPCPort())));
        if (fConnected) break;
        if (fWait)
            MilliSleep(1000);
        else
            throw runtime_error("couldn't connect to server");
    } while (fWait);

    // HTTP basic authentication
    string strUserPass64 = EncodeBase64(mapArgs["-rpcuser"] + ":" + mapArgs["-rpcpassword"]);
    map<string, string> mapRequestHeaders;
    mapRequestHeaders["Authorization"] = string("Basic ") + strUserPass64;

    // Send request
    string strRequest = JSONRPCRequest(strMethod, params, 1);
    string strPost = HTTPPost(strRequest, mapRequestHeaders);
    stream << strPost << std::flush;

    // Receive HTTP reply status
    int nProto = 0;
    int nStatus = ReadHTTPStatus(stream, nProto);

    // Receive HTTP reply message headers and body
    map<string, string> mapHeaders;
    string strReply;
    ReadHTTPMessage(stream, mapHeaders, strReply, nProto);

    if (nStatus == HTTP_UNAUTHORIZED)
        throw runtime_error("incorrect rpcuser or rpcpassword (authorization failed)");
    else if (nStatus >= 400 && nStatus != HTTP_BAD_REQUEST && nStatus != HTTP_NOT_FOUND && nStatus != HTTP_INTERNAL_SERVER_ERROR)
        throw runtime_error(strprintf("server returned HTTP error %d", nStatus));
    else if (strReply.empty())
        throw runtime_error("no response from server");

    // Parse reply
    UniValue valReply(UniValue::VSTR);
    if (!valReply.read(strReply))
        throw runtime_error("couldn't parse reply from server");
    const UniValue& reply = valReply.get_obj();
    if (reply.empty())
        throw runtime_error("expected reply to have result, error and id properties");

    return reply;
}

class CRPCConvertParam
{
public:
    std::string methodName;            // method whose params want conversion
    int paramIdx;                      // 0-based idx of param to convert
};

static const CRPCConvertParam vRPCConvertParams[] =
{
    { "stop", 0 },
    { "getaddednodeinfo", 0 },
    { "sendtoaddress", 1 },
    { "settxfee", 0 },
    { "getstakedbyaddress", 1 },
    { "getreceivedbyaddress", 1 },
    { "getreceivedbyaccount", 1 },
    { "listreceivedbyaddress", 0 },
    { "listreceivedbyaddress", 1 },
    { "listreceivedbyaccount", 0 },
    { "listreceivedbyaccount", 1 },
    { "getnotarytransaction", 1 },
    { "getnotarytransaction", 2 },
    { "getbalance", 1 },
    { "getbalance", 2 },
    { "getblock", 1 },
    { "getblock", 2 },
    { "getblockbynumber", 0 },
    { "getblockbynumber", 1 },
    { "getblockbynumber", 2 },
    { "getblockhash", 0 },
    { "move", 2 },
    { "move", 3 },
    { "sendfrom", 2 },
    { "sendfrom", 3 },
    { "listtransactions", 1 },
    { "listtransactions", 2 },
    { "listbalances", 0 },
    { "listbalances", 1 },
    { "listbalances", 2 },
    { "listaccounts", 0 },
    { "walletpassphrase", 1 },
    { "walletpassphrase", 2 },
    { "getblocktemplate", 0 },
    { "listsinceblock", 1 },
    { "sendalert", 2 },
    { "sendalert", 3 },
    { "sendalert", 4 },
    { "sendalert", 5 },
    { "sendalert", 6 },
    { "sendmany", 1 },
    { "sendmany", 2 },
    { "reservebalance", 0 },
    { "reservebalance", 1 },
    { "addmultisigaddress", 0 },
    { "addmultisigaddress", 1 },
    { "listunspent", 0 },
    { "listunspent", 1 },
    { "listunspent", 2 },
    { "listunspent", 3 },
    { "getrawtransaction", 1 },
    { "createrawtransaction", 0 },
    { "createrawtransaction", 1 },
    { "signrawtransaction", 1 },
    { "signrawtransaction", 2 },
    { "keypoolrefill", 0 },
    { "importprivkey", 2 },
    { "importwallet", 2 },
    { "dumpbootstrap", 1 },
    { "dumpbootstrap", 2 },
    { "validateoutputs", 0 },
    { "listclamours", 0 },
    { "listclamours", 1 },
    { "getsupport", 0 },
    { "getsupport", 1 },
    { "getsupport", 2 },
};

class CRPCConvertTable
{
private:
    std::set<std::pair<std::string, int> > members;

public:
    CRPCConvertTable();

    bool convert(const std::string& method, int idx) {
        return (members.count(std::make_pair(method, idx)) > 0);
    }
};

CRPCConvertTable::CRPCConvertTable()
{
    const unsigned int n_elem =
        (sizeof(vRPCConvertParams) / sizeof(vRPCConvertParams[0]));

    for (unsigned int i = 0; i < n_elem; i++) {
        members.insert(std::make_pair(vRPCConvertParams[i].methodName,
                                      vRPCConvertParams[i].paramIdx));
    }
}

static CRPCConvertTable rpcCvtTable;

/** Non-RFC4627 JSON parser, accepts internal values (such as numbers, true, false, null)
 * as well as objects and arrays.
 */
UniValue ParseNonRFCJSONValue(const std::string& strVal)
{
    UniValue jVal;
    if (!jVal.read(std::string("[")+strVal+std::string("]")) ||
        !jVal.isArray() || jVal.size()!=1)
        throw runtime_error(string("Error parsing JSON:")+strVal);
    return jVal[0];
}

// Convert strings to command-specific RPC representation
UniValue RPCConvertValues(const std::string &strMethod, const std::vector<std::string> &strParams)
{
    UniValue params(UniValue::VARR);

    for (unsigned int idx = 0; idx < strParams.size(); idx++) {
        const std::string& strVal = strParams[idx];

        // insert string value directly
        if (!rpcCvtTable.convert(strMethod, idx)) {
            params.push_back(strVal);
        } else {
            // parse string as JSON, insert bool/number/object/etc. value
            params.push_back(ParseNonRFCJSONValue(strVal));
        }

    }

    return params;
}

int CommandLineRPC(int argc, char *argv[])
{
    string strPrint;
    int nRet = 0;
    try
    {
        // Skip switches
        while (argc > 1 && IsSwitchChar(argv[1][0]))
        {
            argc--;
            argv++;
        }

        // Method
        if (argc < 2)
            throw runtime_error("too few parameters");
        string strMethod = argv[1];

        // Parameters default to strings
        std::vector<std::string> strParams(&argv[2], &argv[argc]);
        UniValue params = RPCConvertValues(strMethod, strParams);

        // Execute
        const UniValue reply = CallRPC(strMethod, params);

        // Parse reply
        const UniValue& result = find_value(reply, "result");
        const UniValue& error  = find_value(reply, "error");

       if (!error.isNull()) 
        {
            // Error
            int code = error["code"].get_int();
            strPrint = "error: " + error.write();
            nRet = abs(code);
            if (error.isObject())
            {
                UniValue errCode = find_value(error, "code");
                UniValue errMsg  = find_value(error, "message");
                strPrint = errCode.isNull() ? "" : "error code: "+errCode.getValStr()+"\n";

                if (errMsg.isStr())
                    strPrint += "error message: "+errMsg.get_str()+"\n";
            }
        }
        else
        {
            // Result
            if (result.isNull())
                strPrint = "";
            else if (result.isStr())
                strPrint = result.get_str();
            else
                strPrint = result.write(2);
        }
    }
    catch (boost::thread_interrupted) {
        throw;
    }
    catch (std::exception& e) {
        strPrint = string("error: ") + e.what();
        nRet = 87;
    }
    catch (...) {
        PrintException(NULL, "CommandLineRPC()");
    }

    if (strPrint != "")
    {
        fprintf((nRet == 0 ? stdout : stderr), "%s\n", strPrint.c_str());
    }
    return nRet;
}
