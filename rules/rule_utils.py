def collect_shape_dirs(config, rules):
    STRIP = "0123456789"
    shapes = set(shape.rstrip(STRIP) for layer in config["layers"].values() for shape in layer.values())
    inputs = {}
    for shape in shapes:
        if shape in ["nuts", "gadm", "lau"]:
            inputs[f"SHAPEINPUT_{shape}"] = getattr(rules, f"administrative_borders_{shape}").output[0]
        else:
            inputs[f"SHAPEINPUT_{shape}"] = config["data-sources"][shape]
    return inputs
