#Execute with your default shell
#i.e.
# . /path/to/script
#      OR (if you want to modify the current busybox .config)
# . /path/to/script bbx
ARCH=arm
#Modify the below variable to your appropriate toolchain path
TOOLC=$HOME/toolchains/arm-uClibc
PATH=$TOOLC/bin:"$PATH"
LD_LIBRARY_PATH=$TOOLC/lib:"$LD_LIBRARY_PATH"
CROSS_COMPILE=arm-linux-
SYSROOT=$TOOLC/arm-rorschack-linux-uclibcgnueabihf/sysroot
CFLAGS="-Wno-error -Os -I$TOOLC/include"
LDFLAGS="-static"
export TOOLC ARCH PATH LD_LIBRARY_PATH CROSS_COMPILE SYSROOT CFLAGS LDFLAGS
export CC="arm-linux-gcc --sysroot=$SYSROOT"
export CPP="arm-linux-cpp --sysroot=$SYSROOT"
export CXX="arm-linux-g++ --sysroot=$SYSROOT"
export LD="arm-linux-ld --sysroot=$SYSROOT"
export AS="arm-linux-as"
export AR="arm-linux-ar"
export RANLIB="arm-linux-ranlib"
if [ _$1 = "_bbx" ];then
	if [ -e $(pwd)/.config ]; then
		sed -i.bak -e "s|.*CONFIG_SYSROOT.*|CONFIG_SYSROOT=\"$SYSROOT\"|" \
		-e "s|.*CONFIG_EXTRA_CFLAGS.*|CONFIG_EXTRA_CFLAGS=\"$CFLAGS\"|" \
		-e "s|.*CONFIG_EXTRA_LDFLAGS.*|CONFIG_EXTRA_LDFLAGS=\"$LDFLAGS\"|" \
		-e "s|.*CONFIG_CROSS_COMPILER_PREFIX.*|CONFIG_CROSS_COMPILER_PREFIX=\"$CROSS_COMPILE\"|" ./.config
	else
		echo "No .config found in current directory!"
	fi
fi
