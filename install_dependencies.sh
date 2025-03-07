# INSTALL DEPENDENCIES
syrahDir=$(pwd)

# ================= GET OS ===========================

if [[ $OSTYPE == "darwin"* ]]; then
  os="mac"
elif [[ $OSTYPE == "linux"* ]]; then
  os="linux"
else
  echo "OS not compatible with installer. Try installing manually or use a Linux or MacOS machine."
  exit 1
fi

# ================= samtools =========================

curl -O -L "https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2"
tar -xjf samtools-1.20.tar.bz2
cd samtools-1.20
./configure --prefix="${syrahDir}/samtools-1.20"
make
make install
export PATH="${syrahDir}/samtools-1.20/bin:$PATH"
cd $syrahDir


# ================= STAR aligner =====================

curl -O -L "https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz"
mv 2.7.10b.tar.gz STAR-2.7.10b.tar.gz
tar -xzf STAR-2.7.10b.tar.gz
cd STAR-2.7.10b/

if [[ $OSTYPE == "darwin"* ]]; then
  # get homebrew installed
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  echo >> /Users/cb2350/.zprofile
  echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/cb2350/.zprofile
  eval "$(/opt/homebrew/bin/brew shellenv)"

  # then gcc and freetype
  brew install gcc
  brew install freetype
  
  # then STAR
  cd source
  make STARforMacStatic CXX=/usr/local/Cellar/gcc/8.2.0/bin/g++-8
  
elif [[ $OSTYPE == "linux"* ]]; then
  cd source
  make STAR
fi

cd $syrahDir

# ================= UMI-tools ========================

curl -O -L "https://github.com/CGATOxford/UMI-tools/archive/refs/tags/1.1.5.tar.gz"
mv 1.1.5.tar.gz UMI-tools-1.1.5.tar.gz
tar -xzf UMI-tools-1.1.5.tar.gz
cd UMI-tools-1.1.5
pip3 install umi_tools==1.1.5

# ================= subread ==========================

if [[ $OSTYPE == "darwin"* ]]; then
  curl -O -L "https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-macOS-x86_64.tar.gz"
  tar -zxf subread-2.0.2-macOS-x86_64.tar.gz
elif [[ $OSTYPE == "linux"* ]]; then
  curl -O -L "https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz"
  tar -zxf subread-2.0.2-Linux-x86_64.tar.gz
fi

echo "ALL DONE!"
