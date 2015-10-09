// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2013 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef _BITCOINRPC_CLIENT_H_
#define _BITCOINRPC_CLIENT_H_ 1

#include <univalue.h>

int CommandLineRPC(int argc, char *argv[]);

UniValue RPCConvertValues(const std::string &strMethod, const std::vector<std::string> &strParams);

#endif
