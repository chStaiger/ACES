#! /bin/bash

set -e

for version in ORIGINAL DETERMINISTIC OPTIONAL_NORMALIZE ONE_THREAD_AND_OPTIONAL_NORMALIZE ; do

    rm -rf PinnacleZ-svn-checkout-2011-04-19
    tar zxf PinnacleZ-svn-checkout-2011-04-19.tar.gz

    # Apply the patch that making PinnacleZ deterministic
    if [ $version = DETERMINISTIC ] ; then
        for PATCH in ForceSingleThreaded.patch InitializeRandomSeeds.patch ; do
            patch -p0 < patches/$PATCH
        done
    fi

    # Apply the patch that making PinnacleZ deterministic
    if [ $version = OPTIONAL_NORMALIZE ] ; then
        for PATCH in AddOptionForNormalizingInput.patch ; do
            patch -p0 < patches/$PATCH
        done
    fi
    
    #added by christine, please check!!!!!!!!!!!!!!!!!!!!!!!!
    # PinnacleZ without z-normalisation and with only one thread
    if [ $version = ONE_THREAD_AND_OPTIONAL_NORMALIZE ]; then
        for PATCH in AddOptionForNormalizingInput.patch ForceSingleThreaded.patch ; do
            patch -p0 < patches/$PATCH
        done
    fi

    cd PinnacleZ-svn-checkout-2011-04-19
    ant
    cd ..

    cp PinnacleZ-svn-checkout-2011-04-19/pinnaclez.jar pinnaclez-$version.jar

    rm -rf PinnacleZ-svn-checkout-2011-04-19

done
