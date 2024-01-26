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

## Examples

### Raw signal acquisition
- Analog voltage signals are converted to digital signals by a National Instruments DAQ
### Behaviour preprocessing
![behaviour preprocessing example](\exampleImages\Behaviour_ex.png "behaviour preprocessing steps")
- Behaviour data is collected at 60Hz.
- Perform light filtering, unwrapping, unit conversion, and resampling as necessary to reduce noise but minimize signal distortion. 
### Ephy preprocessing 
- Electrophysiological data is collected at 5kHz, resample to 1kHz for usability and lightly filter noise without distorting biological signals ocurring at a ms resolution, detect action potentials (spikes) and convert to instantaneous firing rate. 
### Example analyses


## Credits
