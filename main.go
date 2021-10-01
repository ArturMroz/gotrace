package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
)

func main() {
	imageHeight := 256
	imageWidth := 256

	file, err := os.Create("./temp.ppm")
	if err != nil {
		log.Fatal(err)
	}

	w := bufio.NewWriter(file)
	fmt.Fprintf(w, "P3\n%d %d\n255\n", imageHeight, imageWidth)

	for i := 0; i < imageHeight; i++ {
		// log.Println("Scanlines remaining:", imageHeight-i)
		for j := 0; j < imageWidth; j++ {
			c := colour{
				float64(i) / float64(imageHeight),
				float64(j) / float64(imageWidth),
				0.3,
			}

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

func (v vec3) Sub(u vec3) vec3 {
	return vec3{
		x: v.x - u.x,
		y: v.y - u.y,
		z: v.z - u.z,
	}
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
