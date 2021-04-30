import pytest

from renewablepotentialslib.rule_utils import collect_shape_dirs


class TestCustomShapes:
    @pytest.fixture
    def rules(self):
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

    @pytest.mark.parametrize("number", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 100, 1000, 5555))
    def test_strip(self, rules, number):
        config = {"layers": {"national": {"foo": "nuts{}".format(number)}}}
        collected = collect_shape_dirs(config, rules)
        assert collected == {"nuts": "you_got_nuts"}

    @pytest.mark.parametrize("number", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 100, 1000, 5555))
    def test_no_strip(self, rules, number):
        config = {"layers": {"national": {"foo": "{0}bar{0}".format(number)}}, "data-sources": {f"{number}bar": "foobar"}}
        collected = collect_shape_dirs(config, rules)
        assert collected == {f"{number}bar": "foobar"}

    def test_standard_shapes(self, rules):
        config = {
            "layers": {
                "national": {"foo": "nuts2", "bar": "nuts3", "baz": "gadm0"},
                "custom": {"foo": "nuts2", "bar": "lau2", "baz": "gadm0"}
            },
            "data-sources": {"my-custom-shape": "data/dir.geojson"}
        }
        collected = collect_shape_dirs(config, rules)
        assert set(collected.keys()) == set(["lau", "nuts", "gadm"])
        for key in ["lau", "nuts", "gadm"]:
            assert collected[key] == f"you_got_{key}"

    def test_new_shapes(self, rules):
        config = {
            "layers": {
                "national": {"foo": "my-custom-shape1", "bar": "nuts3", "baz": "my-custom-shape2"},
                "custom": {"foo": "your-custom-shape2", "bar": "my-custom-shape12", "baz": "gadm0"}
            },
            "data-sources": {"my-custom-shape": "mydata/dir.geojson", "your-custom-shape": "yourdata/dir.geojson"}
        }
        collected = collect_shape_dirs(config, rules)
        assert set(collected.keys()) == set(["my-custom-shape", "your-custom-shape", "nuts", "gadm"])
        for key in ["my", "your"]:
            assert collected[f"{key}-custom-shape"] == f"{key}data/dir.geojson"
        for key in ["nuts", "gadm"]:
            assert collected[key] == f"you_got_{key}"

    def test_will_fail(self, rules):
        config = {
            "layers": {
                "national": {"foo": "my-custom-shape", "bar": "nuts3"},
            },
            "data-sources": {"your-custom-shape": "yourdata/dir.geojson"}
        }
        with pytest.raises(KeyError) as excinfo:
            collected = set(collect_shape_dirs(config, rules).keys())

        assert "my-custom-shape" in str(excinfo.value)
