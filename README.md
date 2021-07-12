# Current Source Density Analysis

(c) 2015-2021 W.P.Karel Zapfe Zaldivar
México.

The present code is distributed under GNU licence v03.
This code is presented as is, no guarantee whatsoever.

This code presents some tools to obtain the CSD from 2D recordings of
Electrophysiological Activity. The code was developed with 
High Density Microelectrode Arrays (MEAs) in mind. 

There are two main strategies, meant for different use cases.
The first one is the common finite diference approach, for 
signals that have been recorded with small errors or errors outside
the region of interest. The advantages of this approach are speed and
conceptual simplicity. The main disadvantage is that amplifies noise
and errors, so it requires denoising mechanisms in order to give
acceptable results.

The second approach is an inverse problem approach, named
kernel current source density analysis, or kCSDA (Potworowski, bla bla).
This second approach is meant for dealing for recordings that have
many failing electrodes near or in the region of interest. It is also
quite usefull for non rectangular grid matrices. This approach does not
need denoising mechanism, as, by construction, ignores much of the noise.
Its main disadvantage is that is conceptually more sofisticated and 
it requires implementation in parallel in order to compute the
necesary oparators in acceptable times.





