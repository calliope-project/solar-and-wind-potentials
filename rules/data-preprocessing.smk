"""This is a Snakemake file defining rules to retrieve raw data from online sources."""
import pycountry

URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"
URL_NUTS = "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_01M_{}_4326.shp.zip"
URL_LAU = "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/COMM-01M-2013-SH.zip"
URL_DEGURBA = "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/DGURBA_2014_SH.zip"
URL_LAND_COVER = "http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip"
URL_PROTECTED_AREAS = "https://www.protectedplanet.net/downloads/WDPA_Feb2019?type=shapefile"
URL_CGIAR_TILE = "http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/"
URL_GMTED_TILE = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Global_tiles_GMTED/075darcsec/mea/"
URL_GADM = "https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/"
URL_BATHYMETRIC = "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/georeferenced_tiff/ETOPO1_Bed_g_geotiff.zip"
URL_POP = "http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GPW4_GLOBE_R2015A/GHS_POP_GPW42015_GLOBE_R2015A_54009_250/V1-0/GHS_POP_GPW42015_GLOBE_R2015A_54009_250_v1_0.zip"


RAW_SETTLEMENT_DATA = "data/esm-100m-2017/ESM_class{esm_class}_100m.tif"
RAW_EEZ_DATA = "data/World_EEZ_v10_20180221/eez_v10.shp"

RESOLUTION_STUDY = (1 / 3600) * 10 # 10 arcseconds
RESOLUTION_SLOPE = (1 / 3600) * 3 # 3 arcseconds

SRTM_X_MIN = 34 # x coordinate of CGIAR tile raster
SRTM_X_MAX = 44
SRTM_Y_MIN = 1 # y coordinate of CGIAR tile raster
SRTM_Y_MAX = 8
GMTED_Y = ["50N", "70N"]
GMTED_X = ["030W", "000E", "030E"]

localrules: raw_gadm_administrative_borders_zipped, raw_protected_areas_zipped,
    raw_nuts_units_zipped, raw_lau_units_zipped, raw_land_cover_zipped,
    raw_land_cover, raw_protected_areas, raw_srtm_elevation_tile_zipped, raw_gmted_elevation_tile,
    raw_bathymetry_zipped, raw_bathymetry, raw_gadm_administrative_borders


rule raw_gadm_administrative_borders_zipped:
    message: "Download administrative borders for {wildcards.country_code} as zip."
    output: protected("data/automatic/raw-gadm/{country_code}.zip")
    shell: "curl -sLo {output} '{URL_GADM}/gadm36_{wildcards.country_code}_gpkg.zip'"


rule raw_gadm_administrative_borders:
    message: "Unzip administrative borders of {wildcards.country_code} as zip."
    input: "data/automatic/raw-gadm/{country_code}.zip"
    output: temp("data/automatic/raw-gadm/gadm36_{country_code}.gpkg")
    shell: "unzip -o {input} -d data/automatic/raw-gadm"


rule administrative_borders_gadm:
    message: "Merge administrative borders of all countries up to layer {params.max_layer_depth}."
    input:
        "src/gadm.py",
        ["data/automatic/raw-gadm/gadm36_{}.gpkg".format(country_code)
            for country_code in [pycountry.countries.lookup(country).alpha_3
                                 for country in config['scope']['countries']]
         ]
    params: max_layer_depth = 3
    output: "build/administrative-borders-gadm.gpkg"
    conda: "../envs/default.yaml"
    shell:
        PYTHON + " {input} {params.max_layer_depth} {output} {CONFIG_FILE}"


rule raw_nuts_units_zipped:
    message: "Download units as zip."
    output:
        protected("data/automatic/raw-nuts-units.zip")
    shell:
        "curl -sLo {output} '{URL_NUTS}'"


rule administrative_borders_nuts:
    message: "Normalise NUTS administrative borders."
    input:
        src = "src/nuts.py",
        zip = rules.raw_nuts_units_zipped.output
    output:
        "build/administrative-borders-nuts.gpkg"
    shadow: "full"
    params:
        year = config['parameters']['nuts-year']
    conda: "../envs/default.yaml"
    shell:
        """
        unzip {input.zip} -d ./build/NUTS_RG_01M_{params.year}_4326/
        {PYTHON} {input.src} to_multipolygon \
        ./build/NUTS_RG_01M_{params.year}_4326 ./build/raw-nuts.gpkg
        {PYTHON} {input.src} normalise ./build/raw-nuts.gpkg {output} {CONFIG_FILE}
        """


rule raw_lau_units_zipped:
    message: "Download LAU units as zip."
    output:
        protected("data/automatic/raw-lau-units.zip")
    shell:
        "curl -sLo {output} '{URL_LAU}'"


rule administrative_borders_lau:
    message: "Normalise LAU administrative borders."
    input:
        src = "src/lau.py",
        zip = rules.raw_lau_units_zipped.output
    output:
        "build/administrative-borders-lau.geojson"
    shadow: "full"
    conda: "../envs/default.yaml"
    shell:
        """
        unzip {input.zip} -d ./build
        {PYTHON} {input.src} merge ./build/COMM_01M_2013_SH/data/COMM_RG_01M_2013.shp \
        ./build/COMM_01M_2013_SH/data/COMM_AT_2013.dbf ./build/raw-lau.gpkg
        {PYTHON} {input.src} identify ./build/raw-lau.gpkg ./build/raw-lau-identified.gpkg
        {PYTHON} {input.src} normalise ./build/raw-lau-identified.gpkg {output} {CONFIG_FILE}
        """


rule raw_land_cover_zipped:
    message: "Download land cover data as zip."
    output: protected("data/automatic/raw-globcover2009.zip")
    shell: "curl -sLo {output} '{URL_LAND_COVER}'"


rule raw_land_cover:
    message: "Extract land cover data as zip."
    input: rules.raw_land_cover_zipped.output
    output: temp("build/GLOBCOVER_L4_200901_200912_V2.3.tif")
    shadow: "minimal"
    shell: "unzip {input} -d ./build/"


rule raw_protected_areas_zipped:
    message: "Download protected areas data as zip."
    output: protected("data/automatic/raw-wdpa.zip")
    shell: "curl -sLo {output} -H 'Referer: {URL_PROTECTED_AREAS}' {URL_PROTECTED_AREAS}"


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
    output:
        protected("data/automatic/raw-srtm/srtm_{x}_{y}.zip")
    shell:
        """
        curl -sLo {output} '{URL_CGIAR_TILE}/srtm_{wildcards.x}_{wildcards.y}.zip'
        """


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
            "base_url": URL_GMTED_TILE,
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
    shell: "curl -sLo {output} '{URL_BATHYMETRIC}'"


rule raw_bathymetry:
    message: "Extract bathymetric data from zip."
    input: rules.raw_bathymetry_zipped.output
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
        "src/estimate_protected_shapes.py",
        rules.raw_protected_areas.output.points
    output:
        temp("build/protected-areas-points-as-circles.geojson")
    conda: "../envs/default.yaml"
    shell:
        PYTHON_SCRIPT_WITH_CONFIG


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
        class50 = RAW_SETTLEMENT_DATA.format(esm_class="50"),
        class40 = RAW_SETTLEMENT_DATA.format(esm_class="40"),
        class41 = RAW_SETTLEMENT_DATA.format(esm_class="41"),
        class45 = RAW_SETTLEMENT_DATA.format(esm_class="45"),
        class30 = RAW_SETTLEMENT_DATA.format(esm_class="30"),
        class35 = RAW_SETTLEMENT_DATA.format(esm_class="35"),
        reference = rules.land_cover_in_europe.output
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
    shell:
        """
        rio warp {input.bathymetry} -o {output} --like {input.reference} --resampling min
        """


rule eez_in_europe:
    message: "Clip exclusive economic zones to study area."
    input: RAW_EEZ_DATA
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
