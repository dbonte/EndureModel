# EndureModel
Documented code of the mechanistic model developed in the 2Seas Endure project. Supplementary material for manuscript 'Biomorphogenic feedbacks and the spatial organisation of a dominant grass steer dune development' by Dries Bonte<sup>1</sup>, Femke Batsleer<sup>1</sup>, Sebastian Dan<sup>2</sup>, Jasmijn Hillaert<sup>1,3</sup>, Hans Matheve<sup>1</sup>, Sam Provoost<sup>3</sup>, Pieter Rauwoens<sup>4</sup>, Valérie Reijers<sup>5</sup>, Glenn Strypsteen<sup>4</sup>, Suzuki Tomohiro<sup>2</sup>, Martijn Vandegehuchte<sup>1</sup>, Ruben Van de Walle<sup>1</sup>, Toon Verwaest<sup>2</sup>


1.	Ghent University, Dept. Biology, K.L. Ledeganckstraat 35, B-9000 Gent, Belgium
2.	Flanders Hydraulics Research, Berchemlei 115, 2140 Antwerp, Belgium
3.	Institute for Nature and Forest research, Herman Teirlinckgebouw, Havenlaan 88 bus 73, 1000 Brussel
4.	KU Leuven, Campus Brugge, Department of Civil Engineering, Spoorwegstraat 12 8200 Brugge
5.	Department of Physical Geography, Faculty of Geosciences, Utrecht University, Utrecht, 3508 TC the Netherlands 


## A detailed overview of the simulation model.

| Parameter | Definition | Unit | Default value |
| :---: | :--- | :---: | :---: |
| a | Sprouting effect | - | 5 | |
| α | Conversion factor from free-wind velocity to shear-wind velocity | - | 0.058 |
| C | Empirical constant to account for the grain size distribution width | - | 1.8 |
| d<sub>n</sub> | Nominal grain size | µm | 335 |
| D<sub>n</sub> | Reference grain size | µm | 250 |
| d<sub>sand</sub> | Bulk density of sand | kg/m² | 1500 |
| g | Gravitational constant | m/s² | 9.81 |
| Г | Roughness factor of vegetation | - | 16 |
| h | Local vegetation height | m | |
| H | Maximum vegetation height | m | 0.9 |
| l<sub>sat</sub> | Saturation length | m | |
| q | Sand flux | kg/m/s | |
| q<sub>s</sub> | Saturated sand flux | kg/m/s | |
| ρ<sub>a</sub> | Air density | kg/m³ | 1.25 |
| ρ<sub>veg</sub> | Local vegetation density | % | |
| s | Side length of one cell | m | 0.2 |
| τ<sup>*</sup> | Maximum shear stress | Pa | |
| τ<sub>s</sub> | Surface shear stress | Pa | |
| u<sup>*</sup> | Shear velocity | m/s | |
| u<sup>*</sup><sub>t</sub> | Shear velocity threshold | m/s | 3.87 * α |
| u<sub>z</sub> | Wind speed at height z | m/s | 6.41 |
| z | Height above the bed at which the wind speed is measured | m | 10 |

Spatial and temporal dimensions
The landscape represents a square grid with each cell having a dimension of 0.20 x 0.20 m². One time step corresponds to one day.
## Sand dynamics
### Wind direction and boundary conditions

Four different wind directions are defined in the model, each corresponding to a side of the landscape. The distribution of wind direction is assigned at the start of a simulation. Per time step, the wind direction is randomly drawn, based on this distribution. The amount of sand, blown into the system from the sea (N), is expressed as a relative percentage of $q_s$. For instance, if sand influx is defined as 0.5. Then, $q$ equals $0.5 \cdot q_s$ when wind blows from the direction of the sea. Southern winds (land) have a small influx of 1% of $q_s$ sand per cell. Lateral winds have an influx which corresponds with the most recent outflux of a lateral wind, so simulating equal incoming as outcoming lateral fluxes. This amount is constantly updated during a simulation. Wind speed is drawn daily from a normal distribution, based on average wind speed and its standard deviation of that month.

Determine shear velocity based on wind velocity (Hoonhout, 2016): 

$$u^* = \alpha \cdot u_z \tag{eq. 1}$$

Determine maximum (unperturbed) shear stress based on formula for a flatbed (Durán et al., 2010):

$$\tau^* = (u^*)^2 \cdot \rho_a \tag{eq. 2}$$

Calculate fraction of shear stress acting on the sand, based on density of local vegetation (Duran and Moore, 2013):

$$\tau_s = \frac{\tau^*}{1 + \Gamma \cdot \rho_{veg}} \tag{eq. 3}$$

Including Venturi effect: per continuous row of marram grass, perpendicular to the wind direction, the total amount of wind shear stress reduction by vegetation is calculated. A fraction (25%) of this total amount is then added to the wind shear stress of the two adjacent cells of this row. 

Recalculate local wind shear velocity based on formula for an unperturbed shear velocity on a flatbed (Durán et al., 2010):

$$u^* = \sqrt{\tau_s/\rho_a} \tag{eq. 4}$$

Define saturated sand flux per location based on Bagnold formula (Bagnold, 1937; Hoonhout, 2016):

$$q_s = C \cdot \frac{\rho_a}{g} \sqrt{d_n/D_n} \cdot (u^* - u^*_t)^3 \tag{eq. 6}$$

Erosion is modelled based on the following formula (Kroy et al., 2002):

$$\Delta q_{erosion} = \frac{1}{l_{sat}} \cdot q \cdot (1 - q/q_s) \tag{eq. 6}$$

In this formula, $l_{sat}$  is assumed to decrease with wind shear velocity by:

$$l_{sat} = 5/u^* \tag{eq. 7}$$

The maximum amount of sand that can be eroded, is the amount of sand present in a cell.

Deposition is modelled according to:

$$\Delta q_{deposition} = 0.5 \cdot (q - q_s) \tag{eq. 8}$$

### Gravity

Maximum angle of repose is 34° (Durán et al., 2010). Each time step, avalanches are simulated in case this angle is exceeded. Then, the excess amount of sand is displaced to one of the neighbouring cells in the direction with the steepest slope. The maximum angle of repose is set to 34° when vegetation is absent (Durán et al., 2010), but increases with vegetation density. As such, avalanches are less prevalent when plant density is high.


### Shelter effects
Based on vegetation height and sand availability, average slopes (along the wind direction) are determined within the landscapes. In case a lee slope is steeper than 14° (Kroy et al., 2002), a new imaginary slope of 14° is drawn. The area which is covered by this new slope is sheltered. Within this sheltered area, no erosion is allowed. 

### Storm events
A storm event might occur in the middle of a simulation. A cliff erosion simulation can be included with marram grass and sand destroyed in the first 5 m of the landscape, closest to the sea.

### Rain events
Rain events might occur with a chance of 20% (prediction of climate change) from the middle of a simulation onwards.

## Marram grass dynamics

Seasonality in marram grass:
 
### Local growth
Marram grass is only able to grow from April to August (during 153 days), according to:

$$\Delta\rho_{veg} = \rho_{veg} \cdot r \cdot (1 - \rho_{veg}/100) \tag{eq. 9}$$

$r$ represents the growth speed of marram grass and depends on the netto amount of sand deposited or removed (by wind and avalanches) during one time step ($\Delta q_{netto}$) (based on (Nolet et al., 2018). In case no netto deposition of sand occured, growthspeed depends on the number of consecutive days without sand deposition ($t_{no\ deposition}$) or too much sand deposition ($t_{too\ much\ deposition}$). To add randomness to the model, an extra value, drawn from a normal distribution with mean 0 and standard deviation 0.001, is added to $r$ per cell per time step.

$$r = \begin{cases} 
          -462.08 \cdot (\Delta q_{netto} - 0.5/152)^2 + 0.005 & \Delta q_{netto} > 0 \\
          -0.05 & t_{no\ deposition} > 10 \\
          -0.02 & t_{too\ much\ deposition} > 10 
       \end{cases}
       + \mathcal{N}(0,0.001^2) \tag{eq. 10}$$

### Lateral growth
During the growth season, marram grass can also grow laterally. The chance of lateral growth ($\gamma$) depends on $\Delta q_{netto}$ by (based on (Nolet et al., 2018):

$$\gamma = -14440 \cdot (x - 0.4/152)^2 + 0.1 + \mathcal{N}(0, 0.005^2) \tag{eq. 11}$$

In case lateral growth is successful, one of the eight neighbouring cells is randomly selected as the direction of lateral growth. If marram density in that cell is below 90%, the percentage of marram grass is increased by 1. To add randomness to the model, an extra value, drawn from a normal distribution with mean 0 and standard deviation 0.005, is added to $\gamma$ per cell per time step.

### Submersion from September to March

The height of the vegetation in a cell is estimated by (Van Westen, 2018):

$$h_1 = \sqrt{\rho_{veg}} \cdot H \tag{eq. 12}$$

During winter and autumn, vegetation height per cell is updated daily based on the amount of sand deposited or eroded ($\Delta q_{netto}$).

$$h_2 = h_1 - \Delta \frac{q_{netto}}{(d_{sand}/s^2)} \tag{eq. 13}$$

Afterwards, local vegetation cover is estimated based on $h_2$ by (Van Westen, 2018):

$$\rho_{veg} = (\frac{h_2}{H})^2 \tag{eq. 14}$$

### Sprout event at the start of spring
Local density of marram grass after sprouting event is determined by following equation:

$$\rho_{(veg,sprout)} = \rho_{veg} \cdot (1 - \Delta h_{winter}) \cdot a \tag{eq. 15}$$

With $\Delta h$ representing netto change in sand height during autumn and winter. Moreover, in case more than 1 m of sand was locally deposited, marram density becomes 0. In case $\Delta h$ is negative, marram density is unchanged. a determines strength of sprouting effect and was set to five based on field observations. 

## References	
	
* Bagnold, A. R. A. 1937. The Transport of Sand by Wind. - The geographical journal 89: 409–438.
* Durán, O. and Moore, L. J. 2013. Vegetation controls on the maximum size of coastal dunes. - Proc. Natl. Acad. Sci. U. S. A. 110: 17217–22.
* Durán, O. et al. 2010. A continuous model for sand dunes: Review, new developments and application to barchan dunes and barchan dune fields. - Earth Surface Processes and Landforms 35: 1591–1600.
* Hesse, P. P. and Simpson, R. L. 2006. Variable vegetation cover and episodic sand movement on longitudinal desert sand dunes. - Geomorphology 81: 276–291.
* Hoonhout, B. and de Vries, S. 2016. A process-based model for aeolian sediment transport and spatiotemporal varying sediment availability. - Journal of Geophysical Research: Earth Surface 121: 1555–1575.
* Kroy, K. et al. 2002. Minimal model for sand dunes. - Physical Review Letters 88: 543011–543014.
* Nolet, C. et al. 2018. UAV-imaging to model growth response of marram grass to sand burial: Implications for coastal dune development. - Aeolian Research 31: 50–61.
* Van Westen, B. 2018. Numerical modelling of aeolian coastal landform development.
* Van Westen, B. et al. 2019. Aeolian Modelling of Coastal Landform Development. - Coastal Sediments 2019: 1354–1364.
