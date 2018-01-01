Clam 1.4.4
=====================

Copyright (c) 2014 Clam Developers


Setup
---------------------
[Clam Core](http://clamclient.com/download) is the original clam client and it builds the backbone of the network. However, it downloads and stores the entire history of Clam transactions (which is currently several GBs); depending on the speed of your computer and network connection, the synchronization process can take anywhere from a few hours to a day or more. Thankfully you only have to do this once. If you would like the process to go faster you can [download the blockchain directly](https://bitcointalk.org/index.php?topic=623147.msg9772191#msg9772191).

Running
---------------------
The following are some helpful notes on how to run Clam on your native platform. 

### Unix

You need the Qt4 run-time libraries to run Clam-Qt. On Debian or Ubuntu:

	sudo apt-get install libqtgui4

Unpack the files into a directory and run:

- bin/32/clam-qt (GUI, 32-bit) or bin/32/clamd (headless, 32-bit)
- bin/64/clam-qt (GUI, 64-bit) or bin/64/clamd (headless, 64-bit)



### Windows

Unpack the files into a directory, and then run clam-qt.exe.

### OSX

Drag Clam-Qt to your applications folder, and then run Clam-Qt.

### Need Help?

for help and more information.
* Ask for help on [#clams](http://webchat.freenode.net?channels=clams) on Freenode. If you don't have an IRC client use [webchat here](http://webchat.freenode.net?channels=clams).
* Ask for help on the [BitcoinTalk](https://bitcointalk.org/) forums, in the [Technical Support board](https://bitcointalk.org/index.php?topic=623147.0).

Building
---------------------
The following are developer notes on how to build Clam on your native platform. They are not complete guides, but include notes on the necessary libraries, compile flags, etc.

- [OSX Build Notes](build-osx.md)
- [Unix Build Notes](build-unix.md)

Development
---------------------
The Clam repo's [root README](https://github.com/nochowderforyou/clams/blob/master/README.md) contains relevant information on the development process and automated testing.

- [Coding Guidelines](coding.md)
- [Release Notes](release-notes.md)
- [Release Process](release-process.md)
- [Translation Process](translation_process.md)


License
---------------------
Distributed under the [MIT software license](http://www.opensource.org/licenses/mit-license.php).
This product includes software developed by the OpenSSL Project for use in the [OpenSSL Toolkit](https://www.openssl.org/). This product includes
cryptographic software written by Eric Young ([eay@cryptsoft.com](mailto:eay@cryptsoft.com)), and UPnP software written by Thomas Bernard.
