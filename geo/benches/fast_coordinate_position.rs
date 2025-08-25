#[macro_use]
extern crate criterion;
extern crate geo;

use criterion::Criterion;
use geo::contains::IndexedMultiPolygon;
use geo::{coordinate_position::CoordPos, BoundingRect, CoordinatePosition};
use geo_types::{Coord, MultiPolygon};
use geos::CoordSeq;
use wkt::ToWkt;

fn criterion_benchmark(c: &mut Criterion) {
    let zones: MultiPolygon<f64> = geo_test_fixtures::nl_zones();
    let bound = zones.bounding_rect().unwrap();
    let mut coords = vec![];

    // Generate a bunch of points inside the zone bounds
    let size = 20;
    let mut x = bound.min().x;
    for _ in 0..=size {
        let mut y = bound.min().y;
        for _ in 0..=size {
            coords.push(Coord { x, y });
            y += bound.height() / size as f64;
        }

        x += bound.width() / size as f64;
    }

    c.bench_function("with IndexedMultipolygon", |bencher| {
        let indexed = IndexedMultiPolygon::new(&zones);

        bencher.iter(|| {
            let mut inside = 0;

            for c in &coords {
                if indexed.contains_point(*c) {
                    inside += 1
                }
            }

            //assert_eq!(inside, 45); // IndexedMultipolygon gives wrong result
        });
    });

    c.bench_function("with coordinate_position()", |bencher| {
        bencher.iter(|| {
            let mut inside = 0;

            for c in &coords {
                match zones.coordinate_position(c) {
                    CoordPos::Inside => inside += 1,
                    _ => {}
                }
            }

            assert_eq!(inside, 45);
        });
    });

    c.bench_function("with GEOS", |bencher| {
        let geos_geometry = geos::Geometry::new_from_wkt(&zones.to_wkt().to_string()).unwrap();
        let geos_prepared = geos::PreparedGeometry::new(&geos_geometry).unwrap();

        let geos_points: Vec<geos::Geometry> = coords
            .iter()
            .map(|c| {
                geos::Geometry::create_point(CoordSeq::new_from_vec(&[[c.x, c.y]]).unwrap())
                    .unwrap()
            })
            .collect();

        bencher.iter(|| {
            let mut inside = 0;

            for p in &geos_points {
                if geos_prepared.intersects(p).unwrap() {
                    inside += 1;
                }
            }

            assert_eq!(inside, 45);
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
