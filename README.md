## Overview

This repository contains the trunk from the previous subversion (svn) repository of the Los Alamos Sea Ice Model, CICE, including release tags through version 5.1.2. 
More recent versions are found in the [CICE](https://github.com/CICE-Consortium/CICE) and [Icepack](https://github.com/CICE-Consortium/Icepack) repositories, which are now maintained by the CICE Consortium.  

## Obtaining an older CICE release

If you expect to make any changes to the code, we recommend that you work in the CICE and Icepack repositories.  Changes made to code in this repository will not be accepted, other than critical bug fixes.

Release tags from svn have been converted to branches in this git repository.  The pull-down list under "Branch: master" on the [CICE-svn-trunk github page](https://github.com/CICE-Consortium/CICE-svn-trunk) shows all of the available options.
Previous CICE releases may be obtained in several different ways:     
1.  using git    
  git clone -b svn/tags/release-5.1.2 https://github.com/CICE-Consortium/CICE-svn-trunk CICE_v5.1.2
2.  using svn    
  svn checkout https://github.com/CICE-Consortium/CICE-svn-trunk/branches/svn/tags/release-5.1.2 cice_v5.1.2   
3.  download a tarball for a particular version    
[how]
4. clone the entire repository using standard git commands, e.g. 
  git clone https://github.com/CICE-Consortium/CICE-svn-trunk

## More information

Documentation is provided in https://github.com/CICE-Consortium/CICE-svn-trunk/cicedoc/cicedoc.pdf, including both a technical description and user guidance.
Earlier versions of this documentation, for each release, is in the [doc/](https://github.com/CICE-Consortium/CICE-svn-trunk/cice/doc/) directory within the CICE code distribution of that release.

The [wiki](https://github.com/CICE-Consortium/CICE-svn-trunk/wiki) pages for each repository contain links to additional information, e.g.    
- more documentation 
- larger files such as the gx1 grid, land mask, and forcing files

The ["About-Us" repository](https://github.com/CICE-Consortium/About-Us) includes background and supporting information about the CICE Consortium, including how to interact with it.    
