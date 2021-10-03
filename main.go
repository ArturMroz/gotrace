package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
)

const (
	// image

	aspectRatio = 16.0 / 9.0
	imageWidth  = 400
	imageHeight = int(imageWidth / aspectRatio)

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
	// lowerLeftCorner  = origin.Sub(horizontal.Div(2)).Sub(vertical.Div(2)).Sub(vec3{0, 0, focalLength})
)

func main() {
	file, err := os.Create("./temp.ppm")
	if err != nil {
		log.Fatal(err)
	}

	w := bufio.NewWriter(file)
	fmt.Fprintf(w, "P3\n%d %d\n255\n", imageWidth, imageHeight)

	for i := imageHeight - 1; i >= 0; i-- {
		// log.Println("Scanlines remaining:", imageHeight-i)
		for j := imageWidth - 1; j >= 0; j-- {
			u := float64(i) / float64(imageHeight)
			v := float64(j) / float64(imageWidth)
			// r := ray{origin, lowerLeftCorner.Add(horizontal.Mulf(v)).Add(vertical.Mulf(u)).Sub(origin)}
			// r := ray{origin, vadd(lowerLeftCorner, horizontal.Mulf(v), vertical.Mulf(u)).Sub(origin)}
			r := ray{origin, lowerLeftCorner.Add2(horizontal.Mulf(v), vertical.Mulf(u)).Sub(origin)}
			c := r.Colour()

			fmt.Fprintln(w, c.String())
		}
	}

	w.Flush()
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

func (r ray) Colour() colour {
	t := hitSphere(point3{0, 0, -1}, 0.5, r)
	if t > 0 {
		n := Sub(r.At(t), vec3{0, 0, -1}.Unit_vector())
		return Mulf(colour{n.x + 1, n.y + 1, n.z + 1}, 0.5)
	}

	unitDirection := r.dir.Unit_vector()
	t = 0.5 * (unitDirection.y + 1.0)
	// c1 := colour{1.0, 1.0, 1.0}.Mulf(1.0 - t)
	// c2 := colour{0.5, 0.7, 1.0}.Mulf(t)
	// return c1.Add(c2)
	//
	c1 := Mulf(colour{1.0, 1.0, 1.0}, 1.0-t)
	c2 := Mulf(colour{0.5, 0.7, 1.0}, t)
	return Add(c1, c2)
}

func hitSphere(center point3, radius float64, r ray) float64 {
	oc := Sub(r.orig, center)
	a := Dot(r.dir, r.dir)
	b := Dot(oc, r.dir) * 2
	c := Dot(oc, oc) - radius*radius
	// oc := r.orig.Sub(center)
	// a := r.dir.Dot(r.dir)
	// b := oc.Dot(r.dir) * 2
	// c := oc.Dot(oc) - radius*radius
	discriminant := b*b - 4*a*c
	if discriminant < 0 {
		return -1.0
	} else {
		return (-b - math.Sqrt(discriminant)) / (2.0 * a)
	}
}
