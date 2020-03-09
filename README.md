[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

<img src="imgs/logo.png" alt="DLCR" height="28" border="0"/><b><br>Efficient detection of Distant Low Complexity Regions</b>

### Installation ###

Cmake is needed for installation (http://www.cmake.org/). The download is directly available from http://www.cmake.org/cmake/ or by
<pre>
sudo apt-get install cmake git
</pre>

For installing DLCR use:
<pre>
git clone https://github.com/cobilab/DLCR.git
cd DLCR/src/
cmake .
make
</pre>

### Running DLCR ###

<pre>
./DLCR -v -t 1 -l 4 -i 10 bacteria.fa
</pre>

### Parameters ###

<pre>
                                                                        
      DLCR: Efficient detection of Distant Low Complexity Regions       
      ===========================================================       
                                                                        
AUTHORS                                                                 
      Tiago Oliveira (tiagomanuel28@gmail.com) and D. Pratas                     
                                                                        
SYNOPSIS                                                                
      ./DLCR [OPTION]... [FILE]                                         
                                                                        
SAMPLE                                                                  
      ./DLCR -v -l 3 -r 512 -t 1 -w 0.025 genome.fa                     
                                                                        
DESCRIPTION                                                             
      Quantification and Localization of                                
      Distant Low Complexity Regions in FASTA files.                    
                                                                        
      -h,  --help                                                       
           usage guide (help menu).                                     
                                                                        
      -V,  --version                                                    
           Display program and version information.                     
                                                                        
      -F,  --force                                                      
           force mode. Overwrites old files.                            
                                                                        
      -v,  --verbose                                                    
           verbose mode (more information).                             
                                                                        
      -t [NUMBER],  --threshold [NUMBER]                                
           Threshold to segment regions (real).                         
                                                                        
      -w [NUMBER],  --weight [NUMBER]                                   
           Weight to use in low-pass filter (real).                     
                                                                        
      -i [NUMBER],  --ignore [NUMBER]                                   
           Ignore lengths of segmented regions below this value.        
                                                                        
      -r [NUMBER],  --region-size [NUMBER]                              
           Region size to ignore while updating the models.             
           The latest region-size is not updated in the models.         
                                                                        
      -p,  --show-parameters                                            
           show parameters of the models for optimization.              
                                                                        
      -s,  --show-levels                                                
           show pre-computed compression levels (configured parameters).
                                                                        
      -l [NUMBER],  --level [NUMBER]                                    
           Compression level (integer).                                 
           Default level: 4.                                           
           It defines compressibility in balance with computational     
           resources (RAM & time). Use -s for levels perception.        
                                                                        
      [FILE]                                                            
           Input sequence filename (to analyze) -- MANDATORY.           
           File to analyze (last argument).                             
                                                                        
COPYRIGHT                                                               
      Copyright (C) 2020, IEETA, University of Aveiro.                  
      This is a Free software, under GPLv3. You may redistribute        
      copies of it under the terms of the GNU - General Public          
      License v3 <http://www.gnu.org/licenses/gpl.html>. There          
      is NOT ANY WARRANTY, to the extent permitted by law. 
</pre>

### Issues ###

For any issue let us know at [issues link](https://github.com/cobilab/DLCR/issues).

### License ###

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

                                                    

