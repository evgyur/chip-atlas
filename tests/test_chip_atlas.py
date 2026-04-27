import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("chip_atlas", ROOT / "bin" / "chip_atlas.py")
chip_atlas = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules["chip_atlas"] = chip_atlas
SPEC.loader.exec_module(chip_atlas)


class ChipAtlasTests(unittest.TestCase):
    def test_parse_tiny_example(self):
        snps = chip_atlas.parse_raw(ROOT / "examples" / "tiny_atlas_example.txt")
        self.assertEqual(snps["rs9939609"], "AA")
        self.assertEqual(snps["rs671"], "GG")

    def test_build_report_contains_plain_sections(self):
        snps = chip_atlas.parse_raw(ROOT / "examples" / "tiny_atlas_example.txt")
        markers = chip_atlas.build_markers(snps)
        html = chip_atlas.render_html(snps, markers)
        self.assertIn("Вес, аппетит", html)
        self.assertIn("Молочные продукты", html)
        self.assertIn("Что сделать сегодня", html)
        self.assertIn("Consumer raw DNA", html)


if __name__ == "__main__":
    unittest.main()
