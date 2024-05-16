#/bin/sh
make
rnalign_HOME=`pwd`
user_home=`echo ~`
envsetnot=`env | awk -F= '{print $1}' | grep "RNAalign_HOME"`
if [ -z "$envsetnot" ];then 
	cat   << EOF >> ${user_home}/.bashrc
export RNAalign_HOME=${rnalign_HOME}
EOF
	source ${user_home}/.bashrc
fi
