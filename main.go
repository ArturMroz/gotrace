package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
)

const (
	// image

	aspectRatio     = 16.0 / 9.0
	imageWidth      = 400
	imageHeight     = int(imageWidth / aspectRatio)
	samplesPerPixel = 50

	// camera

	viewportHeight = 2.0
	viewportWidth  = aspectRatio * viewportHeight
	focalLength    = 1.0
)

var (
	origin          = colour{0, 0, 0}
	horizontal      = vec3{viewportWidth, 0, 0}
	vertical        = vec3{0, viewportHeight, 0}
	lowerLeftCorner = origin.Sub2(horizontal.Div(2), vertical.Div(2), vec3{0, 0, focalLength})
	world           = hittableList{}
)

func main() {
	file, err := os.Create("./render.ppm")
	if err != nil {
		log.Fatal(err)
	}

	w := bufio.NewWriter(file)
	fmt.Fprintf(w, "P3\n%d %d\n255\n", imageWidth, imageHeight)

	world.objects = append(world.objects, sphere{point3{0, -100.5, -1}, 100})
	world.objects = append(world.objects, sphere{point3{0, 0, -1}, 0.5})

	for i := imageHeight - 1; i >= 0; i-- {
		// log.Println("Scanlines remaining:", imageHeight-i)
		for j := imageWidth - 1; j >= 0; j-- {
			c := colour{}
			for s := 0; s <= samplesPerPixel; s++ {
				u := (float64(i) + rand.Float64()) / float64(imageHeight)
				v := (float64(j) + rand.Float64()) / float64(imageWidth)
				// r := ray{origin, lowerLeftCorner.Add(horizontal.Mulf(v)).Add(vertical.Mulf(u)).Sub(origin)}
				// r := ray{origin, vadd(lowerLeftCorner, horizontal.Mulf(v), vertical.Mulf(u)).Sub(origin)}
				r := ray{origin, lowerLeftCorner.Add2(horizontal.Mulf(v), vertical.Mulf(u)).Sub(origin)}
				c = Add(c, RayColour(r, world))
				// c := RayColour(r, sphere{point3{0, 0, -1}, 0.5})

			}
			fmt.Fprintln(w, c.String2(samplesPerPixel))
		}
	}

	w.Flush()
}

func clamp(x, min, max float64) float64 {
	if x < min {
		return min
	}
	if x > max {
		return max
	}
	return x
}

type vec3 struct {
	x, y, z float64
}

func (v vec3) Add(u vec3) vec3 {
	return vec3{
		x: v.x + u.x,
		y: v.y + u.y,
		z: v.z + u.z,
	}
}

func Add(v, u vec3) vec3 {
	return vec3{
		x: v.x + u.x,
		y: v.y + u.y,
		z: v.z + u.z,
	}
}

func vadd(vs ...vec3) vec3 {
	res := vec3{0, 0, 0}
	for _, v := range vs {
		res.x += v.x
		res.y += v.y
		res.z += v.z
	}
	return res
}

func (v *vec3) Addp(u vec3) {
	v.x -= u.x
	v.y -= u.y
	v.z -= v.z
}

func (v vec3) Add2(us ...vec3) vec3 {
	for _, u := range us {
		v.x += u.x
		v.y += u.y
		v.z += u.z
	}
	return v
}

func (v vec3) Sub(u vec3) vec3 {
	return vec3{
		x: v.x - u.x,
		y: v.y - u.y,
		z: v.z - u.z,
	}
}

func Sub(v, u vec3) vec3 {
	return vec3{
		x: v.x - u.x,
		y: v.y - u.y,
		z: v.z - u.z,
	}
}

func (v vec3) Sub2(us ...vec3) vec3 {
	for _, u := range us {
		v.x -= u.x
		v.y -= u.y
		v.z -= u.z
	}
	return v
}

func (v vec3) Mul(u vec3) vec3 {
	return vec3{
		x: v.x * u.x,
		y: v.y * u.y,
		z: v.z * u.z,
	}
}

func Neg(v vec3) vec3 {
	return vec3{
		x: -v.x,
		y: -v.y,
		z: -v.z,
	}
}

func (v vec3) Mulf(t float64) vec3 {
	return vec3{
		x: v.x * t,
		y: v.y * t,
		z: v.z * t,
	}
}

func Mulf(v vec3, t float64) vec3 {
	return vec3{
		x: v.x * t,
		y: v.y * t,
		z: v.z * t,
	}
}

func (v vec3) Div(t float64) vec3 {
	return v.Mulf(1 / t)
}

func Divf(v vec3, t float64) vec3 {
	return v.Mulf(1 / t)
}

func (v vec3) Len() float64 {
	return math.Sqrt(v.Len_squared())
}

func (v vec3) Len_squared() float64 {
	return v.x*v.x + v.y*v.y + v.z*v.z
}

func (v vec3) Dot(u vec3) float64 {
	return v.x*u.x + v.y*u.y + v.z*u.z
}

func Dot(v, u vec3) float64 {
	return v.x*u.x + v.y*u.y + v.z*u.z
}

func (v vec3) Cross(u vec3) vec3 {
	return vec3{
		x: v.y*u.z - v.z*u.y,
		y: v.z*u.x - v.x*u.z,
		z: v.x*u.y - v.y*v.x,
	}
}

func (v vec3) Unit_vector() vec3 {
	return v.Div(v.Len())
}

func (c colour) String() string {
	return fmt.Sprintf("%d %d %d", int(c.x*255.999), int(c.y*255.999), int(c.z*255.999))
}

func (c colour) String2(samplesPerPixel int) string {
	// TODO that's too much logic for a printing func
	r := c.x
	g := c.y
	b := c.z

	scale := 1.0 / float64(samplesPerPixel)
	r *= scale
	g *= scale
	b *= scale

	return fmt.Sprintf("%d %d %d",
		int(256*clamp(r, 0.0, 0.999)),
		int(256*clamp(g, 0.0, 0.999)),
		int(256*clamp(b, 0.0, 0.999)),
	)
}

type point3 = vec3
type colour = vec3

type ray struct {
	orig point3
	dir  vec3
}

func (r ray) At(t float64) point3 {
	// return r.orig.Add(r.dir.Mulf(t))
	return Add(r.orig, Mulf(r.dir, t))
}

func RayColour(r ray, hittable hittable) colour {
	// var rec hitRecord

	// rec := &hitRecord{}
	// if hittable.hit1(r, 0, math.Inf(1), rec) {

	if rec, isHit := hittable.hit2(r, 0, math.Inf(1)); isHit {
		// return colour{1, 0, 0}

		// log.Println("in ray colour", rec)
		return Mulf(Add(rec.normal, colour{1, 1, 1}), .5)
	}
	// if t > 0 {
	// 	n := Sub(r.At(t), vec3{0, 0, -1}.Unit_vector())
	// 	return Mulf(colour{n.x + 1, n.y + 1, n.z + 1}, 0.5)
	// }

	unitDirection := r.dir.Unit_vector()
	t := 0.5 * (unitDirection.y + 1.0)
	// c1 := colour{1.0, 1.0, 1.0}.Mulf(1.0 - t)
	// c2 := colour{0.5, 0.7, 1.0}.Mulf(t)
	// return c1.Add(c2)
	//
	c1 := Mulf(colour{1.0, 1.0, 1.0}, 1.0-t)
	c2 := Mulf(colour{0.5, 0.7, 1.0}, t)
	return Add(c1, c2)
}

// func hitSphere(center point3, radius float64, r ray) float64 {
// 	oc := Sub(r.orig, center)
// 	a := r.dir.Len_squared()
// 	halfB := Dot(oc, r.dir)
// 	c := oc.Len_squared() - radius*radius
// 	discriminant := halfB*halfB - a*c

// 	if discriminant < 0 {
// 		return -1
// 	} else {
// 		return (-halfB - math.Sqrt(discriminant)) / a
// 	}
// }

type hittable interface {
	hit1(r ray, tMin, tMax float64, rec *hitRecord) bool
	hit2(r ray, tMin, tMax float64) (hitRecord, bool)
}

type hitRecord struct {
	p           point3
	normal      vec3
	t           float64
	isFrontFace bool
}

func (h *hitRecord) setFaceNormal(r ray, outwardNormal vec3) {
	h.isFrontFace = Dot(r.dir, outwardNormal) < 0
	if h.isFrontFace {
		h.normal = outwardNormal
	} else {
		h.normal = Neg(outwardNormal)
	}
}

type sphere struct {
	center point3
	radius float64
}

func (s sphere) hit1(r ray, tMin, tMax float64, rec *hitRecord) bool {
	oc := Sub(r.orig, s.center)
	a := r.dir.Len_squared()
	halfB := Dot(oc, r.dir)
	c := oc.Len_squared() - s.radius*s.radius
	discriminant := halfB*halfB - a*c

	if discriminant < 0 {
		return false
	}

	sqrtd := math.Sqrt(discriminant)
	root := (-halfB - sqrtd) / a
	if tMin > root || root > tMax {
		root := (-halfB + sqrtd) / a
		if tMin > root || root > tMax {
			return false
		}
	}

	rec.t = root
	rec.p = r.At(rec.t)
	outwardNormal := Divf(Sub(rec.p, s.center), s.radius)
	rec.setFaceNormal(r, outwardNormal)

	// log.Println(rec)

	return true
}
func (s sphere) hit2(r ray, tMin, tMax float64) (rec hitRecord, isHit bool) {
	oc := Sub(r.orig, s.center)
	a := r.dir.Len_squared()
	halfB := Dot(oc, r.dir)
	c := oc.Len_squared() - s.radius*s.radius
	discriminant := halfB*halfB - a*c

	if discriminant < 0 {
		return rec, false
	}

	sqrtd := math.Sqrt(discriminant)
	root := (-halfB - sqrtd) / a
	if tMin > root || root > tMax {
		root := (-halfB + sqrtd) / a
		if tMin > root || root > tMax {
			return rec, false
		}
	}

	rec.t = root
	rec.p = r.At(rec.t)
	outwardNormal := Divf(Sub(rec.p, s.center), s.radius)
	rec.setFaceNormal(r, outwardNormal)

	// log.Println(rec)

	return rec, true
}

type hittableList struct {
	objects []hittable
}

func (hl hittableList) hit2(r ray, tMin, tMax float64) (rec hitRecord, hitAnything bool) {
	closestSoFar := tMax
	isHit := false
	tempRec := hitRecord{}

	for _, o := range hl.objects {
		tempRec, isHit = o.hit2(r, tMin, closestSoFar)
		if isHit {
			hitAnything = true
			closestSoFar = tempRec.t
			rec = tempRec
		}
	}

	return rec, hitAnything
}

func (hl hittableList) hit1(r ray, tMin, tMax float64, rec *hitRecord) bool {
	// tempRec := &hitRecord{}
	hitAnything := false
	closestSoFar := tMax

	for _, o := range hl.objects {
		if o.hit1(r, tMin, closestSoFar, rec) {
			hitAnything = true
			closestSoFar = rec.t
			// rec = tempRec
			// rec.t = tempRec.t
			// rec.p = tempRec.p
			// rec.normal = tempRec.normal

		}
	}

	return hitAnything
}

func DegreesToRadians(degrees float64) float64 {
	return degrees * math.Pi / 180
}

// func RandomFloat(min, max) float64
