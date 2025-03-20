# INSTALL DEPENDENCIES
syrahDir=$(pwd)
user_name=$(whoami)

# ================= CHECK OS =========================

if [[ $OSTYPE == "darwin"* ]]; then
  os="mac"
elif [[ $OSTYPE == "linux"* ]]; then
  os="linux"
else
  echo "OS not compatible with installer. Try installing manually or use a Linux or MacOS machine."
  exit 1
fi

# ================= DBSCAN ===========================

echo "install.packages('dbscan',repos='http://cran.us.r-project.org')" | R --vanilla

# ================= samtools =========================

curl -O -L "https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2"
tar -xjf samtools-1.20.tar.bz2
cd samtools-1.20
./configure --prefix="${syrahDir}/samtools-1.20"
make
make install

samtools_path="${syrahDir}/samtools-1.20/bin"
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

  # then mas, xcode
  brew install mas
  mas install 497799835

  # then gcc and freetype
  brew install gcc
  brew install freetype
  
  # then STAR
  cd source
  make STARforMacStatic CXX=/usr/local/Cellar/gcc/8.2.0/bin/g++-8
  star_path="${syrahDir}/STAR-2.7.10b/bin/MacOSX_x86_64"
  
elif [[ $OSTYPE == "linux"* ]]; then
  cd source
  make STAR
  star_path="${syrahDir}/STAR-2.7.10b/bin/Linux_x86_64"
fi

cd $syrahDir

# ================= UMI-tools ========================

curl -O -L "https://github.com/CGATOxford/UMI-tools/archive/refs/tags/1.1.5.tar.gz"
mv 1.1.5.tar.gz UMI-tools-1.1.5.tar.gz
tar -xzf UMI-tools-1.1.5.tar.gz
cd UMI-tools-1.1.5
pip3 install umi_tools==1.1.5
umitools_path="/Users/${user_name}/Library/Python/3.9/lib/python/site-packages"
cd $syrahDir


# ================= subread ==========================

if [[ $OSTYPE == "darwin"* ]]; then
  curl -O -L "https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-macOS-x86_64.tar.gz"
  tar -zxf subread-2.0.2-macOS-x86_64.tar.gz
  subread_path="${syrahDir}/subread-2.0.2-macOS-x86_64/bin"
elif [[ $OSTYPE == "linux"* ]]; then
  curl -O -L "https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz"
  tar -zxf subread-2.0.2-Linux-x86_64.tar.gz
  subread_path="${syrahDir}/subread-2.0.2-Linux-x86_64/bin"
fi

# ================= ADD PATHS ========================

cd $syrahDir

echo "$star_path" > paths.txt
echo "$samtools_path" >> paths.txt
echo "$umitools_path" >> paths.txt
echo "$subread_path" >> paths.txt

echo ""
echo "Paths to install directories have been added to paths.txt and will be added"
echo "to PATH when Syrah is initialized. Feel free to add additional paths (one"
echo "per line) to paths.txt as needed."

echo ""
echo "ALL DONE!"
