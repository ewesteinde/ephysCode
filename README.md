# ephysCode
Personal electrophysiology acquisition, preprocessing &amp; analysis code

## Usage

Code base with 3 major functions:
1. To perform system integration necessary to acquire electrophysiology neural and behaviour signals across a wide range of experiment types and equipment.
    - Contained in ExperimentCode and ExperimentTypes folders
3. To preprocess raw digital signals representing neural activity and behavioural variables of interest.
    - Contained in Analysis_code/General
4. To analyze these signals to understand the link between neural activity and behaviour.
    - Contained in Analysis_code
    - Unique to each cell type and experiment, a lot of test code & WIP, see WesteindeWilson2024 for examples.
      
### Raw signal acquisition

Communication between devices mediated by the Data Acquisition Toolbox with National Instruments support.
  
### Behaviour preprocessing

![behaviour preprocessing steps](https://github.com/ewesteinde/ephysCode/blob/main/exampleImages/Behaviour_ex.png "behaviour preprocessing steps")
Behaviour data is collected at 60Hz.
Perform light filtering, unwrapping, unit conversion, and resampling as necessary to reduce noise but minimize signal distortion. 

### Electrophysiology preprocessing

Detect action potentials (spikes)
![ephys preprocessing step 1](https://github.com/ewesteinde/ephysCode/blob/main/exampleImages/spikedetection_5s.png "spike detection")
Convert binary spikes to instantaneous firing rate, example shown uses a gaussian kernal.
![ephys preprocessing step 2](https://github.com/ewesteinde/ephysCode/blob/main/exampleImages/spikekernal_spikes_5s.png "firing rate conversion")
Electrophysiological data is collected at 5kHz, resampled to 1kHz for usability and lightly filtered noise without distorting biological signals then detect action potentials (spikes) and convert to instantaneous firing rate. 
