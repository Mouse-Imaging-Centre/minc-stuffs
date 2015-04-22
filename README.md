minc-stuffs is repository that houses scripts and bits of code that are useful for performing calculations on and otherwise extracting information from minc files. Contributions are welcome from the entire minc community. 

This code is licensed under the 3-clause BSD License. 
http://opensource.org/licenses/BSD-3-Clause

Installing from github:
-----------------------
<pre><code>
git clone --recursive https://github.com/Mouse-Imaging-Centre/minc-stuffs.git
<\pre><\code>
or
<pre><code>
git clone https://github.com/Mouse-Imaging-Centre/minc-stuffs.git
cd minc-stuffs
git submodule update --init --recursive
<\pre><\code>

To build and install the perl and src code:
------------------------------------------
./autogen.sh
./configure
make
make install

To build and install the python scripts:
---------------------------------------
python setup.py install

Installation Notes:
-------------------
1. If the minc2 libraries are not in your standard path, you will need to include --with-build-path=/path/with/minc-libraries when you run configure.
2. As is standard with configure files, the default installation is in /usr/local/. To install elsewhere, specify --prefix=DIR when you run configure.
3. The default python installation is also /usr/local/. To install elsewhere, specify --home=DIR or --prefix=DIR when running setup.py 

In the future, autogen and configure may be updated to automatically run python setup.py install. 

More information about the package can be found here:

https://wiki.phenogenomics.ca/display/MICePub/minc-stuffs+(formerly+mice-minc-tools)
