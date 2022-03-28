# HOOMD HPS slab builder 
## Code repository to setup slab simulations using the Hydropathy Scale (HPS) Coarse Grained model from a one letter amino acid sequence 
### Authors: Roshan M Regy [1], Wenwei Zheng [2], Jeetain Mittal [1] 
### [1] Department of Chemical and Biomolecular Engineering, Lehigh University 
### [2] College of Integrative Sciences and Arts, Arizona State University
---
### Citations 
1. Regy, R. M.; Zheng, W.; Mittal, J. Theory of biological phase separation, Liquid-Liquid Phase Coexistence and Membraneless Organelles. in Liquid-Liquid Phase Coexistence and Membraneless Organelles (ed. Keating, C. D.) (2020).
2. Dignon, G. L., Zheng, W., Kim, Y. C., Best, R. B. & Mittal, J. Sequence determinants of protein phase behavior from a coarse-grained model. PLoS Comput. Biol. 14, e1005941 (2018).
3. Dignon, G. L., Zheng, W., Best, R. B., Kim, Y. C. & Mittal, J. Relation between single-molecule properties and phase behavior of intrinsically disordered proteins. Proc. Natl. Acad. Sci. 115, 9929–9934 (2018).

--- 
### Prerequisites
The codes provided here were used and tested with the open source MD simulation package [HOOMD-Blue v2.9.0](http://glotzerlab.engin.umich.edu/Downloads/hoomd/) and [gsd v2.1.2](https://gsd.readthedocs.io/en/stable/). To run these codes successfully one would require these aforementioned packages installed. To use HPS interactions between amino acids we also use [azplugins](https://github.com/mphowardlab/azplugins/) which is a plugin for HOOMD-Blue. 

---
### Usage 
1. Create folder with Simulation_with_HOOMD_release.py, stats_module.dat and one letter amino acid sequence file (.dat)
2. Move into above created folder and do the following, 
```
python Simulation_with_HOOMD_release.py seq.dat  
```
---
### Examples
Examples for using this code for simulating FUS LC can be found in directory "FUS_LC_example". On running the code it produces the following files, 

1. ``` start.gsd ``` : contains a file with a single chain constructed from the one letter amino acid sequence in the .dat file passed to the script at runtime 
2. ``` resize.gsd ``` : contains a file with the replicated system (100 chains) reduced to smaller box size (15x15x15nm)
3. ``` box2slab_extend.gsd ``` : contains a file with the energy minimized system in a box extended to the slab configuration (15x15x280nm)
4. ``` Production_dump.dcd ``` : DCD trajectory file with the slab simulated over 1 microsecond 
---
### License and disclaimer
Redistribution and use of this software in source and binary forms, with or without modification, are permitted provided that this statement and the following disclaimer are retained.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# slab_builder
