"""Tests for i.atcorr2.

Standalone tests (no dependency on i.atcorr):
  - Output is physically valid reflectance
  - Higher AOD yields lower BOA reflectance for a dark surface
  - Lower visibility (more aerosol) yields lower BOA reflectance
  - Uniform visibility raster equals visibility_val scalar
  - Uniform elevation raster equals target_elevation scalar
  - Sensor/band SRF wavelength auto-lookup succeeds
  - LUT saved to disk and reloaded gives identical output
  - Uniform AOD raster equals aod_val scalar
  - Active GRASS region is not modified after module run
  - Different aerosol models produce different outputs
  - Different atmosphere models produce different outputs

Cross-comparison with i.atcorr (skipped when i.atcorr is absent):
  Both modules are given the same 6S atmospheric conditions (TM geometry,
  us62 atmosphere, continental aerosol, 40 km visibility) and the same
  reflectance input.  Because they use different 6S engine versions
  (6SV original C++ vs 6SV2.1) and i.atcorr2 uses LUT interpolation,
  exact agreement is not expected; the accepted tolerance is 15 % of the
  mean output value.

  C1. Monochromatic scene: both means agree within 15 %
  C2. Aerosol loading direction: both modules show clear-sky > hazy
  C3. Elevation effect direction: both modules show a nonzero change
      of consistent sign when target is raised from 0 to 1 km
  C4. Visibility raster in i.atcorr ≈ visibility_val= in i.atcorr2
"""

import os
import tempfile
import unittest

import grass.script as gs
from grass.gunittest.case import TestCase
from grass.gunittest.main import test


# ── Shared helper ─────────────────────────────────────────────────────────────

def _write_6s_params(path, *, iwave=-1, wavelength_um=0.660,
                     vis_km=40, target_elev_km=0, atmo=6, aerosol=1):
    """Write a minimal i.atcorr-compatible 6S parameter file.

    Geometry: TM (code 7), June 15 10:00 UTC, lon=0 lat=45 N.
    Spectral: monochromatic (iwave=-1) with explicit wavelength by default.

    Parameters
    ----------
    path : str
        Destination file path.
    iwave : int
        Spectral condition code.  Use -1 for monochromatic (wavelength_um
        is then written on the following line).
    wavelength_um : float
        Wavelength in µm (only used when iwave=-1).
    vis_km : float
        Meteorological visibility in km.
    target_elev_km : float
        Target elevation in km above sea level (written as negative per 6S
        convention).
    atmo : int
        Atmospheric model code (6 = US Standard 1962).
    aerosol : int
        Aerosol model code (1 = continental).
    """
    with open(path, "w") as fh:
        fh.write("7\n")                              # TM geometry code
        fh.write("6 15 10.00 0.0 45.0\n")            # month day hh lon lat
        fh.write(f"{atmo}\n")                        # atmospheric model
        fh.write(f"{aerosol}\n")                     # aerosol model
        fh.write(f"{float(vis_km):.2f}\n")           # visibility km
        fh.write(f"{-abs(float(target_elev_km)):.3f}\n")  # elevation (−km)
        fh.write("-1000\n")                          # satellite sensor height
        if iwave == -1:
            fh.write("-1\n")
            fh.write(f"{wavelength_um:.4f}\n")       # monochromatic wavelength
        else:
            fh.write(f"{iwave}\n")


# ── Standalone tests ─────────────────────────────────────────────────────────

class TestIAtcorr2Standalone(TestCase):
    """i.atcorr2 self-consistency and physical-monotonicity tests.

    All tests use a synthetic 10×10 uniform reflectance raster (ρ = 0.1)
    with the -r flag (input already in reflectance units).  A common set
    of options is shared via :meth:`_common_kwargs` and individual tests
    override only the parameters under test.
    """

    REFL_MAP = "test_atcorr2_refl"
    VIS_MAP  = "test_atcorr2_vis"
    ELEV_MAP = "test_atcorr2_elev"
    AOD_MAP  = "test_atcorr2_aod"

    # Approximate SZA for Jun 15, 10:00 UTC, lat=45 N
    SZA = 32.6
    # Red wavelength matching Sentinel-2A B4 for sensor-lookup tests
    WL  = 0.665

    MODULE = "i.atcorr2"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule("g.region", n=10, s=0, e=10, w=0, rows=10, cols=10)
        # Uniform 10 % reflectance
        cls.runModule("r.mapcalc",
                      expression=f"{cls.REFL_MAP} = 0.1",
                      overwrite=True)
        # Uniform visibility raster 40 km
        cls.runModule("r.mapcalc",
                      expression=f"{cls.VIS_MAP} = 40.0",
                      overwrite=True)
        # Uniform elevation raster 500 m
        cls.runModule("r.mapcalc",
                      expression=f"{cls.ELEV_MAP} = 500.0",
                      overwrite=True)
        # Uniform AOD raster 0.2
        cls.runModule("r.mapcalc",
                      expression=f"{cls.AOD_MAP} = 0.2",
                      overwrite=True)

    @classmethod
    def tearDownClass(cls):
        cls.runModule("g.remove", type="raster", flags="f",
                      pattern="test_atcorr2_*,atcorr2_out_*",
                      quiet=True)
        cls.del_temp_region()

    # ── Helper ──────────────────────────────────────────────────────────── #

    def _common_kwargs(self, output, **overrides):
        """Return a keyword-argument dict for a typical i.atcorr2 invocation.

        The returned dict uses Sentinel-2A B4 (red, ~0.665 µm) with scene-
        average AOD=0.1, H2O=2.0 g/cm², US-62 atmosphere, continental
        aerosol, reflectance input and integer range [0, 1].
        """
        kw = dict(
            input=self.REFL_MAP,
            output=output,
            sensor="sentinel2a",
            band="B4",
            sza=self.SZA,
            atmosphere="us62",
            aerosol="continental",
            aod_val=0.1,
            h2o_val=2.0,
            aod="0.0,0.1,0.2,0.5",
            h2o="1.0,2.0,3.5",
            range="0,1",
            flags="r",
            overwrite=True,
        )
        kw.update(overrides)
        return kw

    @staticmethod
    def _mean(map_name):
        return gs.parse_command("r.univar", map=map_name, format="json")["mean"]

    # ── 1. Basic output range ────────────────────────────────────────────── #

    def test_basic_correction_output_range(self):
        """BOA reflectance from a clear-sky correction must be in [0, 1.5]."""
        out = "atcorr2_out_basic"
        self.assertModule(self.MODULE, **self._common_kwargs(out))
        stats = gs.parse_command("r.univar", map=out, format="json")
        self.assertGreaterEqual(stats["min"], 0.0,
                                "BOA reflectance below 0 — unphysical")
        self.assertLessEqual(stats["max"], 1.5,
                             "BOA reflectance > 1.5 — unphysical")

    # ── 2. AOD sensitivity ───────────────────────────────────────────────── #

    def test_aod_higher_gives_lower_boa(self):
        """For ρ_toa=0.1, BOA at AOD=0.40 must be lower than at AOD=0.05."""
        out_lo = "atcorr2_out_aod_lo"
        out_hi = "atcorr2_out_aod_hi"
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_lo, aod_val=0.05))
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_hi, aod_val=0.40))
        self.assertGreater(self._mean(out_lo), self._mean(out_hi),
                           "Higher AOD did not yield lower BOA reflectance")

    # ── 3. Visibility direction ──────────────────────────────────────────── #

    def test_visibility_lower_gives_lower_boa(self):
        """vis=8 km (hazy) must yield lower BOA than vis=60 km (clear)."""
        out_clear = "atcorr2_out_vis_clear"
        out_hazy  = "atcorr2_out_vis_hazy"
        base = dict(
            input=self.REFL_MAP,
            sza=self.SZA,
            atmosphere="us62",
            aerosol="continental",
            sensor="sentinel2a",
            band="B4",
            h2o_val=2.0,
            aod="0.0,0.1,0.2,0.5",
            h2o="1.0,2.0,3.5",
            range="0,1",
            flags="r",
            overwrite=True,
        )
        self.assertModule(self.MODULE,
                          **{**base, "output": out_clear,
                             "visibility_val": 60})
        self.assertModule(self.MODULE,
                          **{**base, "output": out_hazy,
                             "visibility_val": 8})
        self.assertGreater(self._mean(out_clear), self._mean(out_hazy),
                           "Clear sky (60 km) did not give higher BOA than hazy (8 km)")

    # ── 4. Uniform visibility raster == visibility_val scalar ────────────── #

    def test_uniform_visibility_raster_equals_scalar(self):
        """Uniform 40-km visibility raster must give same output as visibility_val=40."""
        out_rast = "atcorr2_out_vis_raster"
        out_scal = "atcorr2_out_vis_scalar"
        base = dict(
            input=self.REFL_MAP,
            sza=self.SZA,
            atmosphere="us62",
            aerosol="continental",
            sensor="sentinel2a",
            band="B4",
            h2o_val=2.0,
            aod="0.0,0.1,0.2,0.5",
            h2o="1.0,2.0,3.5",
            range="0,1",
            flags="r",
            overwrite=True,
        )
        self.assertModule(self.MODULE,
                          **{**base, "output": out_rast,
                             "visibility": self.VIS_MAP})
        self.assertModule(self.MODULE,
                          **{**base, "output": out_scal,
                             "visibility_val": 40})
        self.assertAlmostEqual(self._mean(out_rast), self._mean(out_scal),
                               delta=1e-4,
                               msg="Visibility raster vs scalar diverge unexpectedly")

    # ── 5. Uniform elevation raster == target_elevation scalar ───────────── #

    def test_uniform_elevation_raster_equals_scalar(self):
        """Elevation raster mean=500 m must equal target_elevation=0.5 km."""
        out_rast = "atcorr2_out_elev_raster"
        out_scal = "atcorr2_out_elev_scalar"
        kw = self._common_kwargs("", target_elevation=0.0)
        kw.pop("output")
        self.assertModule(self.MODULE,
                          **{**kw, "output": out_rast,
                             "elevation": self.ELEV_MAP})
        self.assertModule(self.MODULE,
                          **{**kw, "output": out_scal,
                             "target_elevation": 0.5})
        self.assertAlmostEqual(self._mean(out_rast), self._mean(out_scal),
                               delta=1e-4,
                               msg="Elevation raster vs target_elevation= scalar diverge")

    # ── 6. Sensor/band SRF wavelength lookup ────────────────────────────── #

    def test_sensor_band_lookup_sentinel2a_b2(self):
        """Sensor+band lookup for Sentinel-2A B2 (blue ~0.492 µm) must succeed."""
        out = "atcorr2_out_s2a_b2"
        self.assertModule(self.MODULE,
                          input=self.REFL_MAP,
                          output=out,
                          sensor="sentinel2a",
                          band="B2",
                          sza=self.SZA,
                          atmosphere="us62",
                          aerosol="continental",
                          aod_val=0.15,
                          h2o_val=2.0,
                          aod="0.0,0.1,0.2,0.5",
                          h2o="1.0,2.0,3.5",
                          range="0,1",
                          flags="r",
                          overwrite=True)
        stats = gs.parse_command("r.univar", map=out, format="json")
        self.assertGreaterEqual(stats["n"], 100,
                                "Output has fewer valid pixels than expected")

    def test_sensor_band_lookup_landsat8(self):
        """Sensor+band lookup for Landsat 8 B4 (red ~0.655 µm) must succeed."""
        out = "atcorr2_out_l8_b4"
        self.assertModule(self.MODULE,
                          input=self.REFL_MAP,
                          output=out,
                          sensor="landsat8",
                          band="LC08_B4",
                          sza=self.SZA,
                          atmosphere="us62",
                          aerosol="continental",
                          aod_val=0.1,
                          h2o_val=2.0,
                          aod="0.0,0.1,0.2,0.5",
                          h2o="1.0,2.0,3.5",
                          range="0,1",
                          flags="r",
                          overwrite=True)
        stats = gs.parse_command("r.univar", map=out, format="json")
        self.assertGreaterEqual(stats["n"], 100,
                                "Output has fewer valid pixels than expected")

    # ── 7. LUT save and reload gives identical output ────────────────────── #

    def test_lut_save_and_reload_identical(self):
        """LUT saved to disk and reloaded must produce bit-identical output."""
        out_fresh  = "atcorr2_out_lut_fresh"
        out_loaded = "atcorr2_out_lut_loaded"
        with tempfile.NamedTemporaryFile(suffix=".lut", delete=False) as tmp:
            lut_path = tmp.name
        try:
            # First run: compute and save LUT
            self.assertModule(self.MODULE,
                              **self._common_kwargs(out_fresh, lut=lut_path))
            # Second run: load the saved LUT
            self.assertModule(self.MODULE,
                              **self._common_kwargs(out_loaded, lut=lut_path))
            self.assertAlmostEqual(self._mean(out_fresh), self._mean(out_loaded),
                                   delta=1e-8,
                                   msg="Fresh LUT and loaded LUT outputs differ")
        finally:
            if os.path.exists(lut_path):
                os.unlink(lut_path)

    # ── 8. Uniform AOD raster == aod_val scalar ──────────────────────────── #

    def test_uniform_aod_raster_equals_scalar(self):
        """Uniform AOD raster (0.2) must give same output as aod_val=0.2."""
        out_rast = "atcorr2_out_aod_raster"
        out_scal = "atcorr2_out_aod_scalar"
        base = dict(
            input=self.REFL_MAP,
            sza=self.SZA,
            atmosphere="us62",
            aerosol="continental",
            sensor="sentinel2a",
            band="B4",
            h2o_val=2.0,
            aod="0.0,0.1,0.2,0.5",
            h2o="1.0,2.0,3.5",
            range="0,1",
            flags="r",
            overwrite=True,
        )
        self.assertModule(self.MODULE,
                          **{**base, "output": out_rast,
                             "aod_map": self.AOD_MAP})
        self.assertModule(self.MODULE,
                          **{**base, "output": out_scal,
                             "aod_val": 0.2})
        self.assertAlmostEqual(self._mean(out_rast), self._mean(out_scal),
                               delta=1e-4,
                               msg="AOD raster vs aod_val= scalar diverge unexpectedly")

    # ── 9. Region not modified ───────────────────────────────────────────── #

    def test_region_not_modified_after_run(self):
        """The active GRASS region must be unchanged after i.atcorr2 exits.

        The module adjusts its region to the input raster internally via
        gs.use_temp_region(); the working region of the calling session
        must be restored when the module returns.
        """
        # Create a sub-region that is smaller than the input raster
        self.use_temp_region()
        self.runModule("g.region", n=6, s=0, e=6, w=0, rows=6, cols=6)
        region_before = gs.region()
        out = "atcorr2_out_region_check"
        self.assertModule(self.MODULE, **self._common_kwargs(out))
        region_after = gs.region()
        self.del_temp_region()
        self.assertEqual(region_after["rows"], region_before["rows"],
                         "Module altered the number of rows in the active region")
        self.assertEqual(region_after["cols"], region_before["cols"],
                         "Module altered the number of cols in the active region")
        self.assertAlmostEqual(region_after["n"], region_before["n"], delta=1e-6,
                               msg="Module altered the north boundary of the region")

    # ── 10. Aerosol model produces different output ──────────────────────── #

    def test_aerosol_model_changes_output(self):
        """Continental and urban aerosol models must yield different BOA."""
        out_cont  = "atcorr2_out_aerosol_cont"
        out_urban = "atcorr2_out_aerosol_urban"
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_cont,  aerosol="continental"))
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_urban, aerosol="urban"))
        self.assertNotAlmostEqual(self._mean(out_cont), self._mean(out_urban),
                                  delta=1e-4,
                                  msg="Continental and urban aerosol produced identical output")

    # ── 11. Atmosphere model produces different output ───────────────────── #

    def test_atmosphere_model_changes_output(self):
        """Tropical and mid-latitude-summer atmospheres must yield different BOA."""
        out_trop = "atcorr2_out_atmo_tropical"
        out_mids = "atcorr2_out_atmo_midsum"
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_trop,  atmosphere="tropical"))
        self.assertModule(self.MODULE,
                          **self._common_kwargs(out_mids, atmosphere="midsum"))
        self.assertNotAlmostEqual(self._mean(out_trop), self._mean(out_mids),
                                  delta=1e-5,
                                  msg="Tropical and midsum atmosphere produced identical output")


# ── Cross-comparison with i.atcorr ────────────────────────────────────────────

class TestIAtcorr2VsAtcorr(TestCase):
    """Compare i.atcorr2 output against the classic i.atcorr.

    Both modules are given the same 6S atmospheric conditions and
    reflectance input (``-r`` flag).  Agreement to within 15 % of the
    mean output value is required; differences arise from the different 6S
    engine versions and LUT interpolation in i.atcorr2.

    All tests in this class are skipped automatically when i.atcorr is
    not installed.
    """

    REFL_MAP = "test_atcorr_cmp_refl"
    SZA      = 32.6     # Jun 15, 10:00 UTC, lat=45 N
    WL       = 0.660    # monochromatic red wavelength
    TOL_REL  = 0.15     # 15 % relative tolerance for cross-module comparison

    @classmethod
    def setUpClass(cls):
        """Create synthetic input and check i.atcorr availability."""
        if not gs.find_program("i.atcorr", "--help"):
            raise unittest.SkipTest(
                "i.atcorr is not installed; skipping cross-comparison tests"
            )
        cls.use_temp_region()
        cls.runModule("g.region", n=10, s=0, e=10, w=0, rows=10, cols=10)
        cls.runModule("r.mapcalc",
                      expression=f"{cls.REFL_MAP} = 0.1",
                      overwrite=True)
        # Shared 6S parameter file: monochromatic red, clear sky, sea level
        fd, cls.params_file = tempfile.mkstemp(suffix=".params",
                                               prefix="atcorr_cmp_")
        os.close(fd)
        _write_6s_params(cls.params_file,
                         iwave=-1, wavelength_um=cls.WL,
                         vis_km=40, target_elev_km=0,
                         atmo=6, aerosol=1)

    @classmethod
    def tearDownClass(cls):
        cls.runModule("g.remove", type="raster", flags="f",
                      pattern="test_atcorr_cmp_*,atcorr_cmp_*,atcorr2_cmp_*",
                      quiet=True)
        if hasattr(cls, "params_file") and os.path.exists(cls.params_file):
            os.unlink(cls.params_file)
        cls.del_temp_region()

    # ── Helpers ─────────────────────────────────────────────────────────── #

    def _run_atcorr(self, output, params_file, range_str="0,1",
                    vis_map=None, elev_map=None):
        """Run the classic i.atcorr with reflectance input."""
        kw = dict(
            input=self.REFL_MAP,
            output=output,
            parameters=params_file,
            range=range_str,
            flags="r",
            overwrite=True,
        )
        if vis_map:
            kw["visibility"] = vis_map
        if elev_map:
            kw["elevation"] = elev_map
        self.assertModule("i.atcorr", **kw)

    def _run_atcorr2(self, output, wavelength=None,
                     aod_val=None, visibility_val=None,
                     vis_map=None, elev_map=None,
                     target_elevation=0.0, params_file=None):
        """Run i.atcorr2 with reflectance input."""
        kw = dict(
            input=self.REFL_MAP,
            output=output,
            sza=self.SZA,
            atmosphere="us62",
            aerosol="continental",
            h2o_val=2.0,
            aod="0.0,0.1,0.2,0.5",
            h2o="1.0,2.0,3.5",
            range="0,1",
            target_elevation=target_elevation,
            flags="r",
            overwrite=True,
        )
        if wavelength is not None:
            kw["wavelength"] = wavelength
        if aod_val is not None:
            kw["aod_val"] = aod_val
        if visibility_val is not None:
            kw["visibility_val"] = visibility_val
        if vis_map:
            kw["visibility"] = vis_map
        if elev_map:
            kw["elevation"] = elev_map
        if params_file:
            kw["parameters"] = params_file
        self.assertModule("i.atcorr2", **kw)

    @staticmethod
    def _mean(map_name):
        return gs.parse_command("r.univar", map=map_name, format="json")["mean"]

    def _within_tolerance(self, m1, m2, label1="i.atcorr", label2="i.atcorr2"):
        tol = self.TOL_REL * max(abs(m1), abs(m2), 1e-6)
        self.assertAlmostEqual(
            m1, m2, delta=tol,
            msg=(f"{label1} mean={m1:.5f} vs {label2} mean={m2:.5f} "
                 f"differ by more than {self.TOL_REL * 100:.0f} %")
        )

    # ── C1. Scene mean agreement ─────────────────────────────────────────── #

    def test_reflectance_scene_means_close(self):
        """Monochromatic red correction: both modules agree within 15 %."""
        out1 = "atcorr_cmp_classic"
        out2 = "atcorr2_cmp_new"
        self._run_atcorr(out1,  params_file=self.params_file)
        self._run_atcorr2(out2, wavelength=self.WL, visibility_val=40)
        self._within_tolerance(self._mean(out1), self._mean(out2))

    # ── C2. Aerosol loading direction ────────────────────────────────────── #

    def test_aod_effect_direction_consistent(self):
        """Both modules must show clear-sky > hazy BOA for the same surface."""
        fd1, p_clear = tempfile.mkstemp(suffix=".params", prefix="atcorr_vis50_")
        fd2, p_hazy  = tempfile.mkstemp(suffix=".params", prefix="atcorr_vis8_")
        os.close(fd1)
        os.close(fd2)
        try:
            _write_6s_params(p_clear, iwave=-1, wavelength_um=self.WL, vis_km=50)
            _write_6s_params(p_hazy,  iwave=-1, wavelength_um=self.WL, vis_km=8)

            self._run_atcorr("atcorr_cmp_aod_clear", params_file=p_clear)
            self._run_atcorr("atcorr_cmp_aod_hazy",  params_file=p_hazy)
            self._run_atcorr2("atcorr2_cmp_aod_clear",
                               wavelength=self.WL, visibility_val=50)
            self._run_atcorr2("atcorr2_cmp_aod_hazy",
                               wavelength=self.WL, visibility_val=8)

            m1_clear = self._mean("atcorr_cmp_aod_clear")
            m1_hazy  = self._mean("atcorr_cmp_aod_hazy")
            m2_clear = self._mean("atcorr2_cmp_aod_clear")
            m2_hazy  = self._mean("atcorr2_cmp_aod_hazy")

            self.assertGreater(m1_clear, m1_hazy,
                               "i.atcorr: clear sky did not yield higher BOA than hazy")
            self.assertGreater(m2_clear, m2_hazy,
                               "i.atcorr2: clear sky did not yield higher BOA than hazy")

            # The magnitudes of the difference must also be in the same ballpark
            self._within_tolerance(m1_clear, m2_clear,
                                   "i.atcorr(clear)", "i.atcorr2(clear)")
            self._within_tolerance(m1_hazy, m2_hazy,
                                   "i.atcorr(hazy)", "i.atcorr2(hazy)")
        finally:
            for p in (p_clear, p_hazy):
                if os.path.exists(p):
                    os.unlink(p)

    # ── C3. Elevation effect direction ───────────────────────────────────── #

    def test_elevation_effect_direction_consistent(self):
        """Both modules must show a consistent effect when target rises 1 km."""
        fd_lo, p_lo = tempfile.mkstemp(suffix=".params", prefix="atcorr_elev0_")
        fd_hi, p_hi = tempfile.mkstemp(suffix=".params", prefix="atcorr_elev1_")
        os.close(fd_lo)
        os.close(fd_hi)
        try:
            _write_6s_params(p_lo, iwave=-1, wavelength_um=self.WL,
                             vis_km=40, target_elev_km=0)
            _write_6s_params(p_hi, iwave=-1, wavelength_um=self.WL,
                             vis_km=40, target_elev_km=1)

            self._run_atcorr("atcorr_cmp_elev_0km",  params_file=p_lo)
            self._run_atcorr("atcorr_cmp_elev_1km",  params_file=p_hi)
            self._run_atcorr2("atcorr2_cmp_elev_0km",
                               wavelength=self.WL, visibility_val=40,
                               target_elevation=0.0)
            self._run_atcorr2("atcorr2_cmp_elev_1km",
                               wavelength=self.WL, visibility_val=40,
                               target_elevation=1.0)

            d_classic = (self._mean("atcorr_cmp_elev_1km")
                         - self._mean("atcorr_cmp_elev_0km"))
            d_new     = (self._mean("atcorr2_cmp_elev_1km")
                         - self._mean("atcorr2_cmp_elev_0km"))

            # Each module must individually show a nonzero elevation effect
            self.assertNotAlmostEqual(
                self._mean("atcorr_cmp_elev_0km"),
                self._mean("atcorr_cmp_elev_1km"),
                delta=1e-4,
                msg="i.atcorr shows no effect from 1 km elevation change")
            self.assertNotAlmostEqual(
                self._mean("atcorr2_cmp_elev_0km"),
                self._mean("atcorr2_cmp_elev_1km"),
                delta=1e-4,
                msg="i.atcorr2 shows no effect from 1 km elevation change")

            # Both must agree on the sign of the elevation effect
            self.assertEqual(
                d_classic > 0, d_new > 0,
                msg=(f"i.atcorr and i.atcorr2 disagree on the sign of the "
                     f"elevation effect: i.atcorr Δ={d_classic:.5f}, "
                     f"i.atcorr2 Δ={d_new:.5f}")
            )
        finally:
            for p in (p_lo, p_hi):
                if os.path.exists(p):
                    os.unlink(p)

    # ── C4. Visibility raster in i.atcorr ≈ visibility_val in i.atcorr2 ── #

    def test_visibility_raster_consistent_across_modules(self):
        """Both modules handle visibility input; i.atcorr2 raster ≈ scalar."""
        vis_rast = "test_atcorr_cmp_vis25"
        self.runModule("r.mapcalc",
                       expression=f"{vis_rast} = 25.0",
                       overwrite=True)
        try:
            # i.atcorr with a per-pixel visibility raster
            self._run_atcorr("atcorr_cmp_vis_rast",
                              params_file=self.params_file,
                              vis_map=vis_rast)
            # i.atcorr2 with a per-pixel visibility raster
            self._run_atcorr2("atcorr2_cmp_vis_rast",
                               wavelength=self.WL,
                               vis_map=vis_rast)
            # i.atcorr2 with a scalar visibility value
            self._run_atcorr2("atcorr2_cmp_vis_val",
                               wavelength=self.WL,
                               visibility_val=25)

            m_a2_rast = self._mean("atcorr2_cmp_vis_rast")
            m_a2_val  = self._mean("atcorr2_cmp_vis_val")
            m_classic = self._mean("atcorr_cmp_vis_rast")

            # i.atcorr2: uniform raster == scalar
            self.assertAlmostEqual(
                m_a2_rast, m_a2_val, delta=1e-4,
                msg="i.atcorr2: visibility raster vs scalar differ")

            # Cross-module: both within 15 % of each other
            self._within_tolerance(m_classic, m_a2_val,
                                   "i.atcorr(vis raster)", "i.atcorr2(vis scalar)")
        finally:
            self.runModule("g.remove", type="raster", flags="f",
                           name=vis_rast, quiet=True)

    # ── C5. Aerosol model direction: same ranking in both modules ────────── #

    def test_aerosol_model_ranking_consistent(self):
        """Both modules must agree on the relative ranking of continental vs maritime."""
        fd_cont, p_cont = tempfile.mkstemp(suffix=".params",
                                           prefix="atcorr_continental_")
        fd_mar,  p_mar  = tempfile.mkstemp(suffix=".params",
                                           prefix="atcorr_maritime_")
        os.close(fd_cont)
        os.close(fd_mar)
        try:
            _write_6s_params(p_cont, iwave=-1, wavelength_um=self.WL,
                             vis_km=40, aerosol=1)   # continental
            _write_6s_params(p_mar,  iwave=-1, wavelength_um=self.WL,
                             vis_km=40, aerosol=2)   # maritime

            self._run_atcorr("atcorr_cmp_aerosol_cont", params_file=p_cont)
            self._run_atcorr("atcorr_cmp_aerosol_mar",  params_file=p_mar)
            self._run_atcorr2("atcorr2_cmp_aerosol_cont",
                               wavelength=self.WL, visibility_val=40,
                               params_file=p_cont)
            self._run_atcorr2("atcorr2_cmp_aerosol_mar",
                               wavelength=self.WL, visibility_val=40,
                               params_file=p_mar)

            m1_cont = self._mean("atcorr_cmp_aerosol_cont")
            m1_mar  = self._mean("atcorr_cmp_aerosol_mar")
            m2_cont = self._mean("atcorr2_cmp_aerosol_cont")
            m2_mar  = self._mean("atcorr2_cmp_aerosol_mar")

            # Each module should produce different output for the two aerosol models
            self.assertNotAlmostEqual(m1_cont, m1_mar, delta=1e-4,
                                      msg="i.atcorr: continental == maritime (unexpected)")
            self.assertNotAlmostEqual(m2_cont, m2_mar, delta=1e-4,
                                      msg="i.atcorr2: continental == maritime (unexpected)")

            # Both modules must agree on the direction of the difference
            self.assertEqual(
                (m1_cont > m1_mar), (m2_cont > m2_mar),
                msg=(f"i.atcorr and i.atcorr2 disagree on continental vs maritime "
                     f"ranking: i.atcorr cont={m1_cont:.5f} mar={m1_mar:.5f}, "
                     f"i.atcorr2 cont={m2_cont:.5f} mar={m2_mar:.5f}")
            )
        finally:
            for p in (p_cont, p_mar):
                if os.path.exists(p):
                    os.unlink(p)


if __name__ == "__main__":
    test()
