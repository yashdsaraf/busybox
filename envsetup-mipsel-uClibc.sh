#Execute with your default shell
#i.e.
# . /path/to/script
#      OR (if you want to modify the current busybox .config)
# . /path/to/script bbx
ARCH=mips
#Modify the below variable to your appropriate toolchain path
TOOLC=$HOME/toolchains/mipsel-uClibc
PATH=$TOOLC/bin:"$PATH"
LD_LIBRARY_PATH=$TOOLC/lib:"$LD_LIBRARY_PATH"
CROSS_COMPILE=mipsel-linux-
SYSROOT=$TOOLC/mipsel-rorschack-linux-uclibc/sysroot
CFLAGS="-Wno-error -Os -I$TOOLC/include"
LDFLAGS="-static"
export TOOLC ARCH PATH LD_LIBRARY_PATH CROSS_COMPILE SYSROOT CFLAGS LDFLAGS
export CC="mipsel-linux-gcc --sysroot=$SYSROOT"
export CPP="mipsel-linux-cpp --sysroot=$SYSROOT"
export CXX="mipsel-linux-g++ --sysroot=$SYSROOT"
export LD="mipsel-linux-ld --sysroot=$SYSROOT"
export AS="mipsel-linux-as"
export AR="mipsel-linux-ar"
export RANLIB="mipsel-linux-ranlib"
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
