# EndureModel
Documented code of the mechanistic model developed in the 2Seas Endure project. Supplementary material for manuscript 'Biomorphogenic feedbacks and the spatial organisation of a dominant grass steer dune development' by Dries Bonte, Femke Batsleer, Sebastian Dan, Jasmijn Hillaert, Hans Matheve, Sam Provoost, Pieter Rauwoens, Valérie Reijers, Glenn Strypsteen, Suzuki Tomohiro, Martijn Vandegehuchte, Ruben Van de Walle, Toon Verwaest


1.	Ghent University, Dept. Biology, K.L. Ledeganckstraat 35, B-9000 Gent, Belgium
2.	Flanders Hydraulics Research, Berchemlei 115, 2140 Antwerp, Belgium
3.	Institute for Nature and Forest research, Herman Teirlinckgebouw, Havenlaan 88 bus 73, 1000 Brussel
4.	KU Leuven, Campus Brugge, Department of Civil Engineering, Spoorwegstraat 12 8200 Brugge
5.	Department of Physical Geography, Faculty of Geosciences, Utrecht University, Utrecht, 3508 TC the Netherlands 


A detailed overview of the simulation model.


Parameter	Definition	Unit	Default value
h	Local vegetation height	m	 
lsat	Saturation length	m	 
q	Sand flux	kg/m/s	 
qs	Saturated sand flux	kg/m/s	 
u*	Shear velocity	m/s	 
ρveg	Local vegetation density	%	 
τs	Surface shear stress	Pa	  
a	Sprouting effect	-	5
C	Empirical constant to account for the grain size distribution width	-	1.8
dn	Nominal grain size	µm	335
Dn	Reference grain size	µm	250
dsand	Bulk density of sand	kg/m2	1500
g	Gravitational constant	m/s2	9.81
H	Maximum vegetation height	m	0.9 m
s	Side length of one cell	m	0.2
u*th	Shear velocity threshold	m/s	3.87 * α
uz	Wind speed at height z	m/s	6.41
z	Height above the bed at which the wind speed is measured	m	10
α	Conversion factor from free-wind velocity to shear-wind velocity	-	0.058
ρa	Air density	kg/ m3	1.25
Г	Roughness factor of vegetation	-	16

Spatial and temporal dimensions
The landscape represents a square grid with each cell having a dimension of 0.20 x 0.20 m². One time step corresponds to one day.
Sand dynamics
	Wind direction and boundary conditions

Four different wind directions are defined in the model, each corresponding to a side of the landscape. The distribution of wind direction is assigned at the start of a simulation. Per time step, the wind direction is randomly drawn, based on this distribution. The amount of sand, blown into the system from the sea (N), is expressed as a relative percentage of qs. For instance, if sand influx is defined as 0.5. Then, q equals 0.5 qs when wind blows from the direction of the sea. Southern winds (land) have an influx of 0 kg sand per cell. Lateral winds have an influx which corresponds with the most recent outflux of a lateral wind, so simulating equal incoming as outcoming lateral fluxes. This amount is constantly updated during a simulation. Wind speed is drawn daily from a normal distribution, based on average wind speed and its standard deviation of that month.

	Determine shear velocity based on wind velocity (Hoonhout, 2016): 

u_*= α u_z 	(eq. 1)

	Determine maximum (unperturbed) shear stress based on formula for a flatbed (Durán et al., 2010).
τ_*= 〖u_*〗^2* ρ_a 	(eq. 2)

	Calculate fraction of shear stress acting on the sand, based on density of local vegetation (Duran and Moore, 2013)
 τ_s=  τ_*/(1+ Г ρ_veg ) 		(eq. 3)	
	Including Venturi effect: per continuous row of marram grass, perpendicular to the wind direction, the total amount of wind shear stress reduction by vegetation is calculated. A fraction (25%) of this total amount is then added to the wind shear stress of the two adjacent cells of this row. 

	Recalculate local wind shear velocity based on formula for an unperturbed shear velocity on a flatbed (Durán et al., 2010):

u_*= √(τ_s/ρ_a ) 	(eq. 4)
	Define saturated sand flux per location based on Bagnold formula (Bagnold, 1937; Hoonhout, 2016)

q_s=C  ρ_a/g  √(d_n/D_n )  〖(u_*-u_(*t))〗^3 	(eq. 5)
	Erosion is modelled based on the following formula (Kroy et al., 2002):

∆q_erosion=  1/l_sat   q (1-  q/q_s )	 (eq. 6)

In this formula, l_sat  is assumed to decrease with wind shear velocity by:

l_sat=  5/u_*  	(eq. 7)
The maximum amount of sand that can be eroded, is the amount of sand present in a cell.

	Deposition is modelled according to:

∆q_deposition=0.5*(q- q_s) 	(eq. 8)


	Gravity

Maximum angle of repose is 34° (Durán et al., 2010). Each time step, avalanches are simulated in case this angle is exceeded. Then, the excess amount of sand is displaced to one of the neighbouring cells in the direction with the steepest slope. The maximum angle of repose is set to 34° when vegetation is absent (Durán et al., 2010), but increases with vegetation density. As such, avalanches are less prevalent when plant density is high.

	Shelter effects
Based on vegetation height and sand availability, average slopes (along the wind direction) are determined within the landscapes. In case a lee slope is steeper than 14° (Kroy et al., 2002), a new imaginary slope of 14° is drawn. The area which is covered by this new slope is sheltered. Within this sheltered area, no erosion is allowed. 
	A storm event might occur in the middle of a simulation. A cliff erosion simulation can be included with marram grass and sand destroyed in the first 5 m of the landscape, closest to the sea.

	Rain events might occur with a chance of 20% (prediction of climate change) from the middle of a simulation onwards.

Note: (12) and (13) are not implemented in the here reported simulations

Marram grass dynamics

Seasonality in marram grass:
 
	Local growth
Marram grass is only able to grow from April to August (during 153 days), according to: 
∆ρ_veg= ρ_veg  r (1-  ρ_veg/100) 	(eq. 9)
r represents the growth speed of marram grass and depends on the netto amount of sand deposited or removed (by wind and avalanches) during one time step (∆q_netto) (based on (Nolet et al., 2018). In case no netto deposition of sand occured, growthspeed depends on the number of consecutive days without sand deposition (t_(no deposition)) or too much sand deposition (t_(too much deposition)). To add randomness to the model, an extra value, drawn from a normal distribution with mean 0 and standard deviation 0.01, is added to r per cell per time step.

r={█(-462.08 (∆q_netto- 0.5/152)^2+0.005 if ∆q_netto  >0 @-0.05 if t_(no deposition)>10@-0.05 if t_(too much deposition)>10)}+N(0,〖0.001〗^2)  	(eq. 10)
	Lateral growth
During the growth season, marram grass can also grow laterally. The chance of lateral growth (γ) depends on ∆q_netto by (based on (Nolet et al., 2018):
γ= -14440 (x-0.4/152)^2+0.1+N(0,〖005〗^2)	(eq. 11)
In case lateral growth is successful, one of the eight neighbouring cells is randomly selected as the direction of lateral growth. If marram density in that cell is below 90%, the percentage of marram grass is increased by 1. To add randomness to the model, an extra value, drawn from a normal distribution with mean 0 and standard deviation 0.05, is added to γ per cell per time step.

	Submersion from September to March

The height of the vegetation in a cell is estimated by (Van Westen, 2018):
h_1= 〖√(ρ_veg )〗^   H 	(eq. 12).


During winter and autumn, vegetation height per cell is updated daily based on the amount of sand deposited or eroded (∆q_netto).
 
h_ 2=h_1-(∆q_netto)/(d_sand⁄s^2 ) 	(eq. 13).

Afterwards, local vegetation cover is estimated based on h_2  by (Van Westen, 2018): 				 
ρ_veg= 〖(h_2/H)〗^2 	(eq. 14).	

	Sprout event at the start of spring:
Local density of marram grass after sprouting event is determined by following equation:

ρ_(veg,sprout)= ρ_veg  (1- ∆h_winter)∙a   	(eq. 15).

With ∆h representing netto change in sand height during autumn and winter. Moreover, in case more than 1 m of sand was locally deposited, marram density becomes 0. In case ∆h is negative, marram density is unchanged. a determines strength of sprouting effect and was set to five based on field observations. 


The code is available on Github: TO BE DONE

References	
	
Bagnold, A. R. A. 1937. The Transport of Sand by Wind. - The geographical journal 89: 409–438.
Durán, O. and Moore, L. J. 2013a. Vegetation controls on the maximum size of coastal dunes. - Proc. Natl. Acad. Sci. U. S. A. 110: 17217–22.
Duran, O. and Moore, L. J. 2013. Vegetation controls on the maximum size of coastal dunes. - Proceedings of the National Academy of Sciences 110: 17217–17222.
Durán, O. et al. 2010. A continuous model for sand dunes: Review, new developments and application to barchan dunes and barchan dune fields. - Earth Surface Processes and Landforms 35: 1591–1600.
Hesse, P. P. and Simpson, R. L. 2006. Variable vegetation cover and episodic sand movement on longitudinal desert sand dunes. - Geomorphology 81: 276–291.
Hoonhout, B. and de Vries, S. 2016. A process-based model for aeolian sediment transport and spatiotemporal varying sediment availability. - Journal of Geophysical Research: Earth Surface 121: 1555–1575.
Kroy, K. et al. 2002. Minimal model for sand dunes. - Physical Review Letters 88: 543011–543014.
Nolet, C. et al. 2018. UAV-imaging to model growth response of marram grass to sand burial: Implications for coastal dune development. - Aeolian Research 31: 50–61.
Van Westen, B. 2018. Numerical modelling of aeolian coastal landform development.
Van Westen, B. et al. 2019. AEOLIAN MODELLING OF COASTAL LANDFORM DEVELOPMENT. - Coastal Sediments 2019: 1354–1364.
