"""Script to extract attributes of features and write to csv."""
import click
from fiona.fio.helpers import obj_gen


@click.command()
@click.option('--attribute', '-a', multiple=True)
def geojson_to_csv(attribute):
    """Extract attributes of features and write to csv."""
    stdin = click.get_text_stream('stdin')

    attribute_names = list(attribute)
    click.echo(",".join(attribute_names))

    source = obj_gen(stdin)
    for i, obj in enumerate(source):
        features = obj.get('features') or [obj]
        for j, feat in enumerate(features):
            values = [str(feat["properties"][attribute_name]) if attribute_name in feat["properties"] else ""
                      for attribute_name in attribute_names]
            click.echo(",".join(values))


if __name__ == "__main__":
    geojson_to_csv()
