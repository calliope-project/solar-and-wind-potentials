import pytest

from rules.rule_utils import collect_shape_dirs


class TestCustomShapes:
    PREFIX = "SHAPEINPUT_"
    @pytest.fixture
    def rules(self):
        # This emulates a snakemake rule object, for lack of access to one in the tests
        class DotDict(dict):
            """dot.notation access to dictionary attributes"""
            def __getattr__(*args):
                val = dict.__getitem__(*args)
                return DotDict(val) if type(val) is dict else val
            __setattr__ = dict.__setitem__
            __delattr__ = dict.__delitem__
        return DotDict({
            "administrative_borders_nuts": {"output": ["you_got_nuts"]},
            "administrative_borders_gadm": {"output": ["you_got_gadm"]},
            "administrative_borders_lau": {"output": ["you_got_lau"]}
        })

    @pytest.fixture
    def standard_shapes(self, rules):
        config = {
            "layers": {
                "national": {"foo": "nuts2", "bar": "nuts3", "baz": "gadm0"},
                "custom": {"foo": "nuts2", "bar": "lau2", "baz": "gadm0"}
            },
            "data-sources": {"my-custom-shape": "data/dir.geojson"}
        }
        return collect_shape_dirs(config, rules)

    @pytest.fixture
    def new_shapes(self, rules):
        config = {
            "layers": {
                "national": {"foo": "my-custom-shape1", "bar": "nuts3", "baz": "my-custom-shape2"},
                "custom": {"foo": "your-custom-shape2", "bar": "my-custom-shape12", "baz": "gadm0"}
            },
            "data-sources": {"my-custom-shape": "mydata/dir.geojson", "your-custom-shape": "yourdata/dir.geojson"}
        }
        return collect_shape_dirs(config, rules)

    @pytest.mark.parametrize("number", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 100, 1000, 5555))
    def test_strip_trailing_number(self, rules, number):
        config = {"layers": {"national": {"foo": "nuts{}".format(number)}}}
        collected = collect_shape_dirs(config, rules)
        assert collected == {f"{self.PREFIX}nuts": "you_got_nuts"}

    @pytest.mark.parametrize("number", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 100, 1000, 5555))
    def test_no_strip(self, rules, number):
        config = {"layers": {"national": {"foo": "{0}bar{0}".format(number)}}, "data-sources": {f"{number}bar": "foobar"}}
        collected = collect_shape_dirs(config, rules)
        assert collected == {f"{self.PREFIX}{number}bar": "foobar"}


    def test_collected_standard_shape_sources(self, standard_shapes):
        assert set(standard_shapes.keys()) == set(f"{self.PREFIX}{i}" for i in ["lau", "nuts", "gadm"])

    def test_collected_standard_shapes(self, standard_shapes):
        for key in ["lau", "nuts", "gadm"]:
            assert standard_shapes[f"{self.PREFIX}{key}"] == f"you_got_{key}"

    def test_collected_new_shape_sources(self, new_shapes):
        assert set(new_shapes.keys()) == set(f"{self.PREFIX}{i}" for i in ["my-custom-shape", "your-custom-shape", "nuts", "gadm"])

    def test_collected_new_shapes_when_defined(self, new_shapes):
        for key in ["my", "your"]:
            assert new_shapes[f"{self.PREFIX}{key}-custom-shape"] == f"{key}data/dir.geojson"

    def test_collected_standard_shapes_when_new_shapes_defined(self, new_shapes):
        for key in ["nuts", "gadm"]:
            assert new_shapes[f"{self.PREFIX}{key}"] == f"you_got_{key}"

    def test_no_valid_data_source_for_shape(self, rules):
        config = {
            "layers": {
                "national": {"foo": "my-custom-shape", "bar": "nuts3"},
            },
            "data-sources": {"your-custom-shape": "yourdata/dir.geojson"}
        }
        with pytest.raises(KeyError) as excinfo:
            collected = set(collect_shape_dirs(config, rules).keys())

        assert "my-custom-shape" in str(excinfo.value)
