# Example configuration

`sample_user_inputs.json` contains the parameter set used for the packaged demonstration run. It includes:

- `Vapp_sweep_V`: `1, 0, -1, -5 V`
- interlayer potential-step sweep: `-2, -0.7, +0.7, +2 V`
- left-layer and interlayer decay-length sweeps: `10, 20, 30 nm`
- diffusion coefficients: `1e-19, 1e-18, 1e-17 m2/s`
- elapsed times: `10 and 60 min`
- the supplied default-example ground-state dipole and polarizability tensor

Load these values in the GUI or copy them into the notebook `CFG` cell. Paths and output directories should be adapted to the local machine.
