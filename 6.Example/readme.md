# 📂 6. Example Dataset (PPPH-UAV)

This directory is intended for storing the example GNSS data and associated products required to run the `proc_example.m` script. 

### ⚠️ Data Missing?
Due to GitHub's file size limits (especially for the `.obs`, `.clk`, and `.atx` files), the raw data and products are hosted on **Zenodo** to ensure permanent and open access.

### 📥 Instructions to Get the Data
1. **Download the dataset:** [https://doi.org/10.5281/zenodo.18981904](https://doi.org/10.5281/zenodo.18981904)
2. **Extract:** Unzip the downloaded `PPPH-UAV_Example_Dataset.zip` file.
3. **Move Files:** Place all extracted files directly into this (`6.Example/`) folder.
4. **Verify:** Ensure the following files are present in this directory:
   - `FLY_0098.obs` (Observation data)
   - `COD0MGXFIN_20223090000_01D_01D_OSB.BIA` (Bias products)
   - `COD0MGXFIN_20223090000_01D_05M_ORB.SP3` (Orbit products)
   - `COD0MGXFIN_20223090000_01D_30S_CLK.CLK` (Clock products)
   - `igs20_2233.atx` (Antenna calibration file)
   - `input_example.csv` (UAV event data)
5. **Run:** Open MATLAB, set the current folder to this directory, and run `proc_example.m`.

---
**Reference:**
For more details, please refer to the main repository [README](../README.md) or the published paper: [10.1007/s10291-026-02051-7](https://doi.org/10.1007/s10291-026-02051-7)
