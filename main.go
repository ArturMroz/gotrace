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
	samplesPerPixel = 5
	maxDepth        = 10
)

func main() {
	file, err := os.Create("./render.ppm")
	if err != nil {
		log.Fatal(err)
	}

	w := bufio.NewWriter(file)
	fmt.Fprintf(w, "P3\n%d %d\n255\n", imageWidth, imageHeight)

	world := hittableList{}

	world.objects = append(world.objects, sphere{point3{0, -100.5, -1}, 100, lambertian{colour{.8, .8, 0}}})

	// world.objects = append(world.objects, sphere{point3{0, 0, -1}, 0.5, lambertian{colour{.7, .3, .3}}})
	world.objects = append(world.objects, sphere{point3{0, 0, -1}, 0.5, dielectric{1.5}})
	world.objects = append(world.objects, sphere{point3{-1, 0, -1}, 0.5, dielectric{1.5}})
	// world.objects = append(world.objects, sphere{point3{-1, 0, -1}, -0.4, dielectric{1.5}})
	// world.objects = append(world.objects, sphere{point3{-1, 0, -1}, 0.5, metal{colour{.8, .8, .8}, .3}})
	world.objects = append(world.objects, sphere{point3{1, 0, -1}, 0.5, metal{colour{.8, .6, .2}, 1}})

	cam := newCamera()

	for i := imageHeight - 1; i >= 0; i-- {
		// log.Println("Scanlines remaining:", imageHeight-i)
		for j := imageWidth - 1; j >= 0; j-- {
			c := colour{}
			for s := 0; s <= samplesPerPixel; s++ {
				u := (float64(j) + rand.Float64()) / float64(imageWidth)
				v := (float64(i) + rand.Float64()) / float64(imageHeight)
				// r := ray{origin, lowerLeftCorner.Add2(horizontal.Mulf(v), vertical.Mulf(u)).Sub(origin)}
				r := cam.getRay(u, v)
				c = Add(c, RayColour(r, world, maxDepth))
			}

			fmt.Fprintln(w, c.WriteColour(samplesPerPixel))
		}
	}

	w.Flush()
	log.Println("Done!")
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

type point3 = vec3
type colour = vec3

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

func Mul(v, u vec3) vec3 {
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
	return math.Sqrt(v.LenSquared())
}

func (v vec3) LenSquared() float64 {
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

func (v vec3) UnitVector() vec3 {
	return v.Div(v.Len())
}

func (v vec3) IsNearZero() bool {
	s := 1e-8
	return math.Abs(v.x) < s && math.Abs(v.y) < s && math.Abs(v.z) < s
}

func RandVec() vec3 {
	return vec3{
		x: rand.Float64(),
		y: rand.Float64(),
		z: rand.Float64(),
	}
}

func RandVecRange(min, max float64) vec3 {
	return vec3{
		x: RandFloat64(-1, 1),
		y: RandFloat64(-1, 1),
		z: RandFloat64(-1, 1),
	}
}

func RandInUnitSphere() vec3 {
	for {
		p := RandVecRange(-1, 1)
		if p.LenSquared() < 1 {
			return p
		}
	}
}

func RandUnitVector() vec3 {
	return RandInUnitSphere().UnitVector()
}

func RandFloat64(min, max float64) float64 {
	return min + (max-min)*rand.Float64()
}

func (c colour) String() string {
	return fmt.Sprintf("%d %d %d", int(c.x*255.999), int(c.y*255.999), int(c.z*255.999))
}

func (c colour) WriteColour(samplesPerPixel int) string {
	r := c.x
	g := c.y
	b := c.z

	// divide the colour by number of samples and gamma correct for gamma=2.0
	scale := 1.0 / float64(samplesPerPixel)
	r = math.Sqrt(r * scale)
	g = math.Sqrt(g * scale)
	b = math.Sqrt(b * scale)

	return fmt.Sprintf("%d %d %d",
		int(256*clamp(r, 0, 0.999)),
		int(256*clamp(g, 0, 0.999)),
		int(256*clamp(b, 0, 0.999)),
	)
}

type ray struct {
	orig point3
	dir  vec3
}

func (r ray) At(t float64) point3 {
	// return r.orig.Add(r.dir.Mulf(t))
	return Add(r.orig, Mulf(r.dir, t))
}

func Reflect(v, n vec3) vec3 {
	// log.Println("ref:", Sub(v, Mulf(n, 2*Dot(v, n))))
	return Sub(v, Mulf(n, 2*Dot(v, n)))
}

func Refract(uv, n vec3, etaiOverEtat float64) vec3 {
	cosTheta := math.Min(Dot(Neg(uv), n), 1.0)
	rOutPerpendicular := Mulf(Add(uv, Mulf(n, cosTheta)), etaiOverEtat)
	rOutParallel := Mulf(n, -math.Sqrt(math.Abs(1.0-rOutPerpendicular.LenSquared())))
	// log.Println(rOutPerpendicular)
	// log.Println(rOutParallel)
	// log.Println("added:", Add(rOutPerpendicular, rOutParallel))
	return Add(rOutPerpendicular, rOutParallel)
}

func Refract2(v vec3, n vec3, nRatio float64) vec3 {
	cosTheta := Dot(Neg(v), n)
	parallel := Mulf(Add(v, Mulf(n, cosTheta)), nRatio)
	sqrt := -math.Sqrt(1.0 - parallel.LenSquared())
	normal := Mulf(n, sqrt)
	return Add(normal, parallel)
}

func RayColour(r ray, hittable hittable, depth int) colour {
	if depth <= 0 {
		return colour{1, 0, 0}
	}

	if rec, isHit := hittable.hit2(r, 0.001, math.Inf(1)); isHit {
		var scattered ray
		var attenuation colour

		if rec.material.Scatter(&r, &rec, &attenuation, &scattered) {
			return Mul(attenuation, RayColour(scattered, world, depth-1))
		}

		return colour{0, 1, 0}
	}

	unitDirection := r.dir.UnitVector()
	t := 0.5 * (unitDirection.y + 1.0)
	c1 := Mulf(colour{1.0, 1.0, 1.0}, 1.0-t)
	c2 := Mulf(colour{0.5, 0.7, 1.0}, t)

	return Add(c1, c2)
}

type hittable interface {
	hit1(r ray, tMin, tMax float64, rec *hitRecord) bool
	hit2(r ray, tMin, tMax float64) (hitRecord, bool)
}

type hitRecord struct {
	p           point3
	normal      vec3
	t           float64
	isFrontFace bool
	material    material
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
	center   point3
	radius   float64
	material material
}

func (s sphere) hit1(r ray, tMin, tMax float64, rec *hitRecord) bool {
	oc := Sub(r.orig, s.center)
	a := r.dir.LenSquared()
	halfB := Dot(oc, r.dir)
	c := oc.LenSquared() - s.radius*s.radius
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

	return true
}

func (s sphere) hit2(r ray, tMin, tMax float64) (rec hitRecord, isHit bool) {
	oc := Sub(r.orig, s.center)
	a := r.dir.LenSquared()
	halfB := Dot(oc, r.dir)
	c := oc.LenSquared() - s.radius*s.radius
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
	rec.material = s.material

	return rec, true
}

type hittableList struct {
	objects []hittable
}

func (hl hittableList) hit2(r ray, tMin, tMax float64) (rec hitRecord, hitAnything bool) {
	closestSoFar := tMax

	for _, o := range hl.objects {
		if tempRec, isHit := o.hit2(r, tMin, closestSoFar); isHit {
			hitAnything = true
			closestSoFar = tempRec.t
			rec = tempRec
		}
	}

	return rec, hitAnything
}

func (hl hittableList) hit1(r ray, tMin, tMax float64, rec *hitRecord) bool {
	hitAnything := false
	closestSoFar := tMax

	for _, o := range hl.objects {
		if o.hit1(r, tMin, closestSoFar, rec) {
			hitAnything = true
			closestSoFar = rec.t
		}
	}

	return hitAnything
}

func DegreesToRadians(degrees float64) float64 {
	return degrees * math.Pi / 180
}

type material interface {
	Scatter(rIn *ray, rec *hitRecord, attenuation *colour, scattered *ray) bool
}

type lambertian struct {
	albedo colour
}

func (l lambertian) Scatter(rIn *ray, rec *hitRecord, attenuation *colour, scattered *ray) bool {
	scatterDirection := Add(rec.normal, RandUnitVector())
	// catch degenerate scatter direction
	if scatterDirection.IsNearZero() {
		scatterDirection = rec.normal
	}

	*scattered = ray{rec.p, scatterDirection}
	*attenuation = l.albedo

	return true
}

type metal struct {
	albedo colour
	fuzz   float64
}

func (m metal) Scatter(rIn *ray, rec *hitRecord, attenuation *colour, scattered *ray) bool {
	reflected := Reflect(rIn.dir.UnitVector(), rec.normal)
	*scattered = ray{rec.p, Add(reflected, Mulf(RandInUnitSphere(), m.fuzz))}
	*attenuation = m.albedo

	return Dot(scattered.dir, rec.normal) > 0
}

type dielectric struct {
	// index of refraction
	ir float64
}

func (d dielectric) Scatter(rIn *ray, rec *hitRecord, attenuation *colour, scattered *ray) bool {
	*attenuation = colour{1, 1, 1}
	var refractionRatio float64
	if rec.isFrontFace {
		refractionRatio = 1.0 / d.ir
	} else {
		refractionRatio = d.ir
	}

	unitDir := rIn.dir.UnitVector()
	cosTheta := math.Min(Dot(Neg(unitDir), rec.normal), 1.0)
	sinTheta := math.Sqrt(1.0 - cosTheta*cosTheta)

	cannotRefract := refractionRatio*sinTheta > 1

	var direction vec3
	if cannotRefract || reflectance(cosTheta, refractionRatio) > rand.Float64() {
		direction = Reflect(unitDir, rec.normal)
	} else {
		direction = Refract2(rIn.dir, rec.normal, refractionRatio)
	}

	// refracted := Refract(rIn.dir.UnitVector(), rec.normal, refractionRatio)
	// refracted := Refract(rIn.dir, rec.normal, refractionRatio)
	// log.Println("refracted in scatter:", refracted)

	*scattered = ray{rec.p, direction}

	return true
}

func reflectance(cosine, refIdx float64) float64 {
	// Shilck's approximation for reflectance
	r0 := (1 - refIdx) / (1 + refIdx)
	r0 *= r0
	return r0 + (1-r0)*math.Pow((1-cosine), 5)
}

type camera struct {
	origin          point3
	lowerLeftCorner point3
	horizontal      vec3
	vertical        vec3
}

func newCamera() camera {
	aspectRatio := 16.0 / 9.0
	viewportHeight := 2.0
	viewportWidth := aspectRatio * viewportHeight
	focalLength := 1.0

	c := camera{
		origin:     colour{0, 0, 0},
		horizontal: vec3{viewportWidth, 0, 0},
		vertical:   vec3{0, viewportHeight, 0},
	}

	c.lowerLeftCorner = c.origin.Sub2(c.horizontal.Div(2), c.vertical.Div(2), vec3{0, 0, focalLength})

	return c
}

func (c camera) getRay(u, v float64) ray {
	dir := c.lowerLeftCorner.Add(c.horizontal.Mulf(u)).Add(c.vertical.Mulf(v)).Sub(c.origin)
	return ray{c.origin, dir}
}
