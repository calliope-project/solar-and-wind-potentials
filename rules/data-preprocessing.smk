"""This is a Snakemake file defining rules to retrieve raw data from online sources."""
import pycountry
from rule_utils import collect_shape_dirs

RESOLUTION_STUDY = (1 / 3600) * 10 # 10 arcseconds
RESOLUTION_SLOPE = (1 / 3600) * 3 # 3 arcseconds

SRTM_X_MIN = 34 # x coordinate of CGIAR tile raster
SRTM_X_MAX = 44
SRTM_Y_MIN = 1 # y coordinate of CGIAR tile raster
SRTM_Y_MAX = 8
GMTED_Y = ["50N", "70N"]
GMTED_X = ["030W", "000E", "030E"]

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"

localrules: raw_gadm_administrative_borders_zipped, raw_protected_areas_zipped,
    raw_lau_units_zipped, raw_land_cover_zipped,
    raw_land_cover, raw_protected_areas, raw_srtm_elevation_tile_zipped, raw_gmted_elevation_tile,
    raw_bathymetry_zipped, raw_bathymetry, raw_gadm_administrative_borders


rule raw_gadm_administrative_borders_zipped:
    message: "Download administrative borders for {wildcards.country_code} as zip."
    params: url = lambda wildcards: config["data-sources"]["gadm"].format(country_code=wildcards.country_code)
    output: protected("data/automatic/raw-gadm/{country_code}.zip")
    shell: "curl -sLo {output} '{params.url}'"


rule raw_gadm_administrative_borders:
    message: "Unzip administrative borders of {wildcards.country_code}."
    input: "data/automatic/raw-gadm/{country_code}.zip"
    output: temp("build/raw-gadm/gadm36_{country_code}.gpkg")
    shell: "unzip -o {input} -d build/raw-gadm"


rule administrative_borders_gadm:
    message: "Merge gadm administrative borders of all countries."
    input:
        ["build/raw-gadm/gadm36_{}.gpkg".format(country_code)
            for country_code in [pycountry.countries.lookup(country).alpha_3
                                 for country in config['scope']['countries']]
         ]
    output: temp("build/raw-gadm/gadm36.gpkg")
    params: crs = config["crs"]
    conda: '../envs/default.yaml'
    shell: "ogrmerge.py -o {output} -f gpkg -src_layer_field_content '{{LAYER_NAME}}' -t_srs {params.crs} -single {input}"


rule administrative_borders_nuts:
    message: "Download NUTS units as GeoJSON."
    output:
        protected("data/automatic/raw-nuts{}-units.geojson".format(config["parameters"]["nuts-year"]))
    params:
        url = config["data-sources"]["nuts"].format(nuts_year=config["parameters"]["nuts-year"])
    shell:
        "curl -sLo {output} '{params.url}'"


rule raw_lau_units_zipped:
    message: "Download LAU units as zip."
    output:
        protected("data/automatic/raw-lau-units.zip")
    params: url = config["data-sources"]["lau"]
    shell:
        "curl -sLo {output} '{params.url}'"


rule raw_lau_units_unzipped:
    message: "Unzip LAU units."
    input:
        zip = rules.raw_lau_units_zipped.output
    output:
        shapes = "build/raw-lau-units/COMM_RG_01M_2013.shp",
        attributes = "build/raw-lau-units/COMM_AT_2013.dbf"
    shell: "unzip -j {input.zip} -d build/raw-lau-units"


rule administrative_borders_lau:
    message: "Normalise LAU administrative borders."
    input:
        src = script_dir + "lau.py",
        shapes = rules.raw_lau_units_unzipped.output.shapes,
        attributes = rules.raw_lau_units_unzipped.output.attributes,
    output:
        temp("build/raw-lau.gpkg")
    shadow: "full"
    conda: "../envs/default.yaml"
    script: "../scripts/lau.py"


rule administrative_borders:
    message: "Normalise all administrative borders."
    input:
        src = script_dir + "administrative_borders.py",
        **collect_shape_dirs(config, rules)
    params:
        crs = config["crs"],
        scope = config["scope"]
    output:
        "build/administrative-borders.gpkg"
    shadow: "full"
    conda: "../envs/default.yaml"
    script: "../scripts/administrative_borders.py"


rule raw_land_cover_zipped:
    message: "Download land cover data as zip."
    output: protected("data/automatic/raw-globcover2009.zip")
    params: url = config["data-sources"]["land_cover"]
    shell: "curl -sLo {output} '{url}'"


rule raw_land_cover:
    message: "Extract land cover data as zip."
    input: rules.raw_land_cover_zipped.output
    output: temp("build/GLOBCOVER_L4_200901_200912_V2.3.tif")
    shadow: "minimal"
    shell: "unzip {input} -d ./build/"


rule raw_protected_areas_zipped:
    message: "Download protected areas data as zip."
    output: protected("data/automatic/raw-wdpa.zip")
    params: url = config["data-sources"]["protected_areas"]
    shell: "curl -sLo {output} -H 'Referer: {params.url}' {params.url}"


rule raw_protected_areas:
    message: "Extract protected areas data as zip."
    input: rules.raw_protected_areas_zipped.output
    output:
        polygons = "build/raw-wdpa-feb2019/WDPA_Feb2019-shapefile-polygons.shp",
        polygon_data = "build/raw-wdpa-feb2019/WDPA_Feb2019-shapefile-polygons.dbf",
        points = "build/raw-wdpa-feb2019/WDPA_Feb2019-shapefile-points.shp"
    shell: "unzip -o {input} -d build/raw-wdpa-feb2019"


rule raw_srtm_elevation_tile_zipped:
    message: "Download SRTM elevation data tile (x={wildcards.x}, y={wildcards.y}) from CGIAR."
    output: protected("data/automatic/raw-srtm/srtm_{x}_{y}.zip")
    params: url = config["data-sources"]["cgiar_tile"]
    shell: "curl -sLo {output} '{params.url}/srtm_{wildcards.x}_{wildcards.y}.zip'"


rule raw_srtm_elevation_tile:
    message: "Unzip SRTM elevation data tile (x={wildcards.x}, y={wildcards.y})."
    input:
        "data/automatic/raw-srtm/srtm_{x}_{y}.zip"
    output:
        temp("build/srtm_{x}_{y}.tif")
    run: # using Python here, because "unzip" issues a warning sometimes causing snakemake to break
        import zipfile
        from pathlib import Path

        file_to_extract = Path(input[0]).with_suffix(".tif").name
        with zipfile.ZipFile(input[0], "r") as zipref:
            zipref.extract(file_to_extract, path="build")


rule raw_srtm_elevation_data:
    message: "Merge all SRTM elevation data tiles."
    input:
        ["build/srtm_{x:02d}_{y:02d}.tif".format(x=x, y=y)
         for x in range(SRTM_X_MIN, SRTM_X_MAX + 1)
         for y in range(SRTM_Y_MIN, SRTM_Y_MAX + 1)
         if not (x is 34 and y in [3, 4, 5, 6])] # these tiles do not exist
    output:
        temp("build/raw-srtm-elevation-data.tif")
    conda: "../envs/default.yaml"
    shell:
        "rio merge {input} {output} --overwrite"


rule raw_gmted_elevation_tile:
    message: "Download GMTED elevation data tile."
    output:
        protected("data/automatic/raw-gmted/raw-gmted-{y}-{x}.tif")
    run:
        url = "{base_url}/{x_inverse}/{y}{x}_20101117_gmted_mea075.tif".format(**{
            "base_url": config["data-sources"]["gmted_tile"],
            "x": wildcards.x,
            "y": wildcards.y,
            "x_inverse": wildcards.x[-1] + wildcards.x[:-1]
        })
        shell("curl -sLo {output} '{url}'".format(**{"url": url, "output": output}))


rule raw_gmted_elevation_data:
    message: "Merge all GMTED elevation data tiles."
    input:
        ["data/automatic/raw-gmted/raw-gmted-{y}-{x}.tif".format(x=x, y=y)
         for x in GMTED_X
         for y in GMTED_Y
         ]
    output:
        temp("build/raw-gmted-elevation-data.tif")
    conda: "../envs/default.yaml"
    shell:
        "rio merge {input} {output} --overwrite"


rule raw_bathymetry_zipped:
    message: "Download bathymetric data as zip."
    output: protected("data/automatic/raw-bathymetric.zip")
    params: url = config["data-sources"]["bathymetric"]
    shell: "curl -sLo {output} '{params.url}'"


rule raw_bathymetry:
    message: "Extract bathymetric data from zip."
    input: rules.raw_bathymetry_zipped.output[0]
    output: temp("build/ETOPO1_Bed_g_geotiff.tif")
    shell: "unzip {input} -d ./build/"


rule elevation_in_europe:
    message: "Merge SRTM and GMTED elevation data and warp/clip to Europe using {threads} threads."
    input:
        gmted = rules.raw_gmted_elevation_data.output,
        srtm = rules.raw_srtm_elevation_data.output
    output:
        temp("build/elevation-europe.tif")
    params:
        srtm_bounds = "{x_min},{y_min},{x_max},60".format(**config["scope"]["bounds"]),
        gmted_bounds = "{x_min},59.5,{x_max},{y_max}".format(**config["scope"]["bounds"])
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    shell:
        """
        rio clip --bounds {params.srtm_bounds} {input.srtm} -o build/tmp-srtm.tif
        rio clip --bounds {params.gmted_bounds} {input.gmted} -o build/tmp-gmted.tif
        rio warp build/tmp-gmted.tif -o build/tmp-gmted2.tif -r {RESOLUTION_SLOPE} \
        --resampling nearest --threads {threads}
        rio merge build/tmp-srtm.tif build/tmp-gmted2.tif {output}
        rm build/tmp-gmted.tif
        rm build/tmp-gmted2.tif
        rm build/tmp-srtm.tif
        """


rule land_cover_in_europe:
    message: "Clip land cover data to Europe."
    input: rules.raw_land_cover.output
    output: "build/land-cover-europe.tif"
    params: bounds = "{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["bounds"])
    conda: "../envs/default.yaml"
    shell: "rio clip {input} {output} --bounds {params.bounds}"


rule slope_in_europe:
    message: "Calculate slope and warp to resolution of study using {threads} threads."
    input:
        elevation = rules.elevation_in_europe.output,
        land_cover = rules.land_cover_in_europe.output
    output:
        "build/slope-europe.tif"
    threads: config["snakemake"]["max-threads"]
    conda: "../envs/default.yaml"
    shell:
        """
        gdaldem slope -s 111120 -compute_edges {input.elevation} build/slope-temp.tif
        rio warp build/slope-temp.tif -o {output} --like {input.land_cover} \
        --resampling max --threads {threads}
        rm build/slope-temp.tif
        """


rule protected_areas_points_to_circles:
    message: "Estimate shape of protected areas available as points only."
    input:
        script_dir + "estimate_protected_shapes.py",
        protected_areas= rules.raw_protected_areas.output.points
    params:
        scope = config["scope"]
    output:
        temp("build/protected-areas-points-as-circles.geojson")
    conda: "../envs/default.yaml"
    script: "../scripts/estimate_protected_shapes.py"


rule protected_areas_in_europe:
    message: "Rasterise protected areas data and clip to Europe."
    input:
        polygons = rules.raw_protected_areas.output.polygons,
        points = rules.protected_areas_points_to_circles.output,
        land_cover = rules.land_cover_in_europe.output
    output:
        "build/protected-areas-europe.tif"
    benchmark:
        "build/rasterisation-benchmark.txt"
    params:
        bounds = "{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["bounds"])
    conda: "../envs/default.yaml"
    shell:
        # The filter is in accordance to the way UNEP-WCMC calculates statistics:
        # https://www.protectedplanet.net/c/calculating-protected-area-coverage
        """
        fio cat --rs --bbox {params.bounds} {input.polygons} {input.points} | \
        fio filter "f.properties.STATUS in ['Designated', 'Inscribed', 'Established'] and \
        f.properties.DESIG_ENG != 'UNESCO-MAB Biosphere Reserve'" | \
        fio collect --record-buffered | \
        rio rasterize --like {input.land_cover} \
        --default-value 255 --all_touched -f "GTiff" --co dtype=uint8 -o {output}
        """


rule settlements:
    message: "Warp settlement data to CRS of study using {threads} threads."
    input:
        class50 = config["data-sources"]["settlement_data"].format(esm_class="50"),
        class40 = config["data-sources"]["settlement_data"].format(esm_class="40"),
        class41 = config["data-sources"]["settlement_data"].format(esm_class="41"),
        class45 = config["data-sources"]["settlement_data"].format(esm_class="45"),
        class30 = config["data-sources"]["settlement_data"].format(esm_class="30"),
        class35 = config["data-sources"]["settlement_data"].format(esm_class="35"),
        reference = rules.land_cover_in_europe.output[0]
    output:
        buildings = "build/esm-class50-buildings.tif",
        urban_greens = "build/esm-class404145-urban-greens.tif",
        built_up = "build/esm-class303550-built-up.tif"
    threads: config["snakemake"]["max-threads"]
    shadow: "minimal"
    conda: "../envs/default.yaml"
    shell:
        """
        rio calc "(+ (+ (read 1) (read 2)) (read 3))" \
        {input.class40} {input.class41} {input.class45} -o build/esm-class404145-temp-not-warped.tif
        rio calc "(+ (+ (read 1) (read 2)) (read 3))" \
        {input.class50} {input.class30} {input.class35} -o build/esm-class303550-temp-not-warped.tif
        rio warp {input.class50} -o {output.buildings} \
        --like {input.reference} --threads {threads} --resampling bilinear
        rio warp build/esm-class404145-temp-not-warped.tif -o {output.urban_greens} \
        --like {input.reference} --threads {threads} --resampling bilinear
        rio warp build/esm-class303550-temp-not-warped.tif -o {output.built_up} \
        --like {input.reference} --threads {threads} --resampling bilinear
        """


rule bathymetry_in_europe:
    message: "Clip bathymetric data to study area and warp to study resolution."
    input:
        bathymetry = rules.raw_bathymetry.output,
        reference = rules.land_cover_in_europe.output
    output:
        "build/bathymetry-in-europe.tif"
    conda: "../envs/default.yaml"
    shell: "rio warp {input.bathymetry} -o {output} --like {input.reference} --resampling min"


rule eez_in_europe:
    message: "Clip exclusive economic zones to study area."
    input: config["data-sources"]["eez_data"]
    output: "build/eez-in-europe.geojson"
    params:
        bounds="{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["bounds"]),
        countries=",".join(["'{}'".format(country) for country in config["scope"]["countries"]])
    conda: "../envs/default.yaml"
    shell:
        """
        fio cat --bbox {params.bounds} {input}\
        | fio filter "f.properties.Territory1 in [{params.countries}]"\
        | fio collect > {output}
        """
