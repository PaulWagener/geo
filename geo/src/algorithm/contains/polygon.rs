use super::{impl_contains_from_relate, impl_contains_geometry_for, Contains};
use crate::geometry::*;
use crate::{GeoFloat, GeoNum};
use crate::{HasDimensions, Relate};
use ordered_float::OrderedFloat;
use std::ops::ControlFlow;

// ┌─────────────────────────────┐
// │ Implementations for Polygon │
// └─────────────────────────────┘
impl<T> Contains<Coord<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, coord: &Coord<T>) -> bool {
        use crate::coordinate_position::{CoordPos, CoordinatePosition};

        self.coordinate_position(coord) == CoordPos::Inside
    }
}

impl<T> Contains<Point<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, p: &Point<T>) -> bool {
        self.contains(&p.0)
    }
}

impl<T> Contains<MultiPoint<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, mp: &MultiPoint<T>) -> bool {
        use crate::coordinate_position::{CoordPos, CoordinatePosition};
        // at least one point must be fully within
        // others can be on the boundary
        mp.iter().any(|p| self.contains(p))
            && mp
                .iter()
                .all(|p| self.coordinate_position(&p.0) != CoordPos::Outside)
    }
}

impl_contains_from_relate!(Polygon<T>, [Line<T>, LineString<T>, Polygon<T>, MultiLineString<T>, MultiPolygon<T>, GeometryCollection<T>, Rect<T>, Triangle<T>]);
impl_contains_geometry_for!(Polygon<T>);

// ┌──────────────────────────────────┐
// │ Implementations for MultiPolygon │
// └──────────────────────────────────┘

impl<T> Contains<Coord<T>> for MultiPolygon<T>
where
    T: GeoNum,
{
    fn contains(&self, coord: &Coord<T>) -> bool {
        self.iter().any(|poly| poly.contains(coord))
    }
}

impl<T> Contains<Point<T>> for MultiPolygon<T>
where
    T: GeoNum,
{
    fn contains(&self, p: &Point<T>) -> bool {
        self.contains(&p.0)
    }
}

impl<T: GeoNum> Contains<MultiPoint<T>> for MultiPolygon<T> {
    fn contains(&self, rhs: &MultiPoint<T>) -> bool {
        if self.is_empty() || rhs.is_empty() {
            return false;
        }
        rhs.iter().all(|point| self.contains(point))
    }
}

impl<F> Contains<Line<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &Line<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<LineString<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &LineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<MultiLineString<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &MultiLineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<Polygon<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &Polygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<MultiPolygon<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &MultiPolygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<GeometryCollection<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &GeometryCollection<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<Rect<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &Rect<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F> Contains<Triangle<F>> for MultiPolygon<F>
where
    F: GeoFloat,
{
    fn contains(&self, rhs: &Triangle<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

// Extension trait for optimized point-in-multipolygon test
pub trait ContainsPointFast<T: GeoFloat + rstar::RTreeNum> {
    /// An optimized point-in-multipolygon test using an interval tree.
    fn contains_point_fast(&self, point: &Point<T>) -> bool;
}

impl<T> ContainsPointFast<T> for MultiPolygon<T>
where
    T: GeoFloat,
{
    fn contains_point_fast(&self, _point: &Point<T>) -> bool {
        todo!()
    }
}

use crate::intersects::value_in_between;
use crate::kernels::{Kernel, Orientation};
use crate::{Coord, LineString, MultiPolygon, Polygon};
use sif_itree::ITree;

struct YIntervalSegment<F: GeoFloat> {
    segment: (Coord<F>, Coord<F>),
    is_exterior: bool,
}

pub struct IndexedMultiPolygon {
    y_interval_tree: ITree<OrderedFloat<f64>, YIntervalSegment<f64>>,
    geometry: MultiPolygon<f64>,
}

impl IndexedMultiPolygon {
    pub fn new(mp: &MultiPolygon) -> Self {
        let mp = mp.clone();
        let mut segments = Vec::new();

        for polygon in &mp.0 {
            Self::add_ring_segments(&mut segments, polygon.exterior(), true);
            for interior in polygon.interiors() {
                Self::add_ring_segments(&mut segments, interior, false);
            }
        }

        Self {
            y_interval_tree: ITree::new(segments),
            geometry: mp,
        }
    }

    fn add_ring_segments(
        segments: &mut Vec<sif_itree::Item<OrderedFloat<f64>, YIntervalSegment<f64>>>,
        ring: &LineString<f64>,
        is_exterior: bool,
    ) {
        for window in ring.coords().collect::<Vec<_>>().windows(2) {
            let (p1, p2) = (window[0], window[1]);
            let y_min = OrderedFloat(p1.y.min(p2.y));
            let y_max = OrderedFloat(p1.y.max(p2.y));
            segments.push((
                y_min..y_max,
                YIntervalSegment {
                    segment: (*p1, *p2),
                    is_exterior,
                },
            ));
        }
    }

    pub fn contains_point(&self, point: Coord<f64>) -> bool {
        // Use winding number algorithm with robust predicates
        // Based on coord_pos_relative_to_ring in coordinate_position.rs
        let mut winding_number = 0;

        let result = self.y_interval_tree.query(
            OrderedFloat(point.y)..OrderedFloat(point.y),
            |(_, segment)| {
                let seg = segment.segment;

                // Apply winding number algorithm using robust predicates
                let mut contribution = 0;
                if seg.0.y <= point.y {
                    if seg.1.y >= point.y {
                        let o = <f64 as GeoNum>::Ker::orient2d(seg.0, seg.1, point);
                        if o == Orientation::CounterClockwise && seg.1.y != point.y {
                            contribution = 1;
                        } else if o == Orientation::Collinear
                            && value_in_between(point.x, seg.0.x, seg.1.x)
                        {
                            // Point on boundary!
                            return ControlFlow::Break(false);
                        }
                    }
                } else if seg.1.y <= point.y {
                    let o = <f64 as GeoNum>::Ker::orient2d(seg.0, seg.1, point);
                    if o == Orientation::Clockwise {
                        contribution = -1;
                    } else if o == Orientation::Collinear
                        && value_in_between(point.x, seg.0.x, seg.1.x)
                    {
                        // Point on boundary!
                        return ControlFlow::Break(false);
                    }
                }

                // Exterior rings contribute positively, interior rings (holes) contribute negatively
                if segment.is_exterior {
                    winding_number += contribution;
                } else {
                    winding_number -= contribution;
                }
                ControlFlow::Continue(())
            },
        );

        match result {
            ControlFlow::Break(r) => r,
            ControlFlow::Continue(()) => winding_number != 0,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{coord, polygon, MultiPoint, Relate};

    fn make_test_pts() -> [Coord<f64>; 7] {
        let pt_a = coord! {x: 0., y: 0.};
        let pt_b = coord! {x: 10., y: 0.};
        let pt_c = coord! {x: 10., y: 10.};
        let pt_d = coord! {x: 0., y: 10.};

        let pt_edge = coord! {x: 0., y: 5.};
        let pt_mid = coord! {x: 5., y: 5.};
        let pt_out = coord! {x: 11., y: 11.};
        [pt_a, pt_b, pt_c, pt_d, pt_edge, pt_mid, pt_out]
    }

    #[test]
    fn test_polygon_should_never_contain_empty_multipoint() {
        let [pt_a, pt_b, pt_c, pt_d, _pt_edge, _pt_mid, _pt_out] = make_test_pts();

        let poly = polygon![pt_a, pt_b, pt_c, pt_d, pt_a];
        let empty: MultiPoint<f64> = MultiPoint::new(Vec::new());

        // contains implementation follows `Relate`` trait
        assert!(!poly.contains(&empty));
        assert!(!poly.relate(&empty).is_contains());
    }

    #[test]
    fn test_polygon_should_contains_multipoint() {
        let [pt_a, pt_b, pt_c, pt_d, pt_edge, pt_mid, _pt_out] = make_test_pts();

        let poly = polygon![pt_a, pt_b, pt_c, pt_d, pt_a];

        // contains requires at least one point fully within the polygon
        let mp_a_mid = MultiPoint::from(vec![pt_a, pt_mid]);
        let mp_bc_mid = MultiPoint::from(vec![pt_a, pt_b, pt_mid]);
        let mp_bc_edge_mid = MultiPoint::from(vec![pt_a, pt_b, pt_edge, pt_mid]);

        assert!(poly.contains(&mp_a_mid));
        assert!(poly.contains(&mp_bc_mid));
        assert!(poly.contains(&mp_bc_edge_mid));
    }

    #[test]
    fn test_polygon_should_not_contains_multipoint_on_edge() {
        let [pt_a, pt_b, pt_c, pt_d, pt_edge, _pt_mid, _pt_out] = make_test_pts();

        let poly = polygon![pt_a, pt_b, pt_c, pt_d, pt_a];

        // contains should return false if all points lie on the boundary of the polygon
        let mp_a = MultiPoint::from(vec![pt_a]);
        let mp_bc = MultiPoint::from(vec![pt_a, pt_b]);
        let mp_bc_edge = MultiPoint::from(vec![pt_a, pt_b, pt_edge]);

        assert!(!poly.contains(&mp_a));
        assert!(!poly.contains(&mp_bc));
        assert!(!poly.contains(&mp_bc_edge));
    }

    #[test]
    fn test_polygon_should_not_contains_out_multipoint() {
        let [pt_a, pt_b, pt_c, pt_d, pt_edge, _pt_mid, pt_out] = make_test_pts();

        let poly = polygon![pt_a, pt_b, pt_c, pt_d, pt_a];

        // contains should return false if any point lies outside the polygon
        let mp_a_out = MultiPoint::from(vec![pt_a, pt_out]);
        let mp_bc_out = MultiPoint::from(vec![pt_a, pt_b, pt_out]);
        let mp_bc_edge_out = MultiPoint::from(vec![pt_a, pt_b, pt_edge, pt_out]);

        assert!(!poly.contains(&mp_a_out));
        assert!(!poly.contains(&mp_bc_out));
        assert!(!poly.contains(&mp_bc_edge_out));
    }
}
