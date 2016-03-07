#! /bin/sh

set -e

BASENAME=`date +"PinnacleZ-svn-checkout-%Y-%m-%d"`
svn co http://chianti.ucsd.edu/svn/csplugins/trunk/ucsd/slotia/pinnaclez "${BASENAME}"
rm -rf `find "${BASENAME}" -type d -name ".svn"`
tar zcvf "${BASENAME}.tar.gz" "${BASENAME}"
rm -rf "${BASENAME}"

BASENAME=`date +"Oiler-svn-checkout-%Y-%m-%d"`
svn co http://chianti.ucsd.edu/svn/csplugins/trunk/ucsd/slotia/oiler "${BASENAME}"
rm -rf `find "${BASENAME}" -type d -name ".svn"`
tar zcvf "${BASENAME}.tar.gz" "${BASENAME}"
rm -rf "${BASENAME}"

BASENAME=`date +"ModLab-svn-checkout-%Y-%m-%d"`
svn co http://chianti.ucsd.edu/svn/csplugins/trunk/ucsd/slotia/modlab "${BASENAME}"
rm -rf `find "${BASENAME}" -type d -name ".svn"`
tar zcvf "${BASENAME}.tar.gz" "${BASENAME}"
rm -rf "${BASENAME}"
