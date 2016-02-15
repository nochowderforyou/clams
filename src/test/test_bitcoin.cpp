
// Copyright (c) 2011-2013 The Bitcoin Core developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#define BOOST_TEST_MODULE Bitcoin Test Suite

#include "db.h"
#include "txdb.h"
#include "main.h"
#include "ui_interface.h"
#include "wallet.h"
#include "util.h"

CWallet* pwalletMain;
CClientUIInterface uiInterface;


#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

extern bool fPrintToConsole;
extern void noui_connect();

struct TestingSetup {
    CCoinsViewDB *pcoinsdbview;
    boost::filesystem::path pathTemp;
    boost::thread_group threadGroup;

    TestingSetup() {
        fPrintToDebugLog = false; // don't want to write to debug.log file
        noui_connect();
        bitdb.MakeMock();
        pathTemp = GetTempPath() / strprintf("test_clam_%lu_%i", (unsigned long)GetTime(), (int)(GetRand(100000)));
        boost::filesystem::create_directories(pathTemp);
        mapArgs["-datadir"] = pathTemp.string();
        pblocktree = new CBlockTreeDB(1 << 20, true);
        pcoinsdbview = new CCoinsViewDB(1 << 23, true);
        pcoinsTip = new CCoinsViewCache(*pcoinsdbview);
        InitBlockIndex();
        bool fFirstRun;
        pwalletMain = new CWallet("wallet.dat");
        pwalletMain->LoadWallet(fFirstRun);
        RegisterWallet(pwalletMain);

    }
    ~TestingSetup()
    {
        threadGroup.interrupt_all();
        threadGroup.join_all();
        delete pwalletMain;
        pwalletMain = NULL;
        delete pcoinsTip;
        delete pcoinsdbview;
        delete pblocktree;
        bitdb.Flush(true);
        boost::filesystem::remove_all(pathTemp);
    }
};

BOOST_GLOBAL_FIXTURE(TestingSetup);
