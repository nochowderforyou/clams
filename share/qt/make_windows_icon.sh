#!/bin/bash
# create multiresolution windows icon
ICON_DST=../../src/qt/res/icons/clam.ico

convert ../../src/qt/res/icons/clam-16.png ../../src/qt/res/icons/clam-32.png ../../src/qt/res/icons/clam-48.png ${ICON_DST}
