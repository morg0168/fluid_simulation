
// INCOMPRESSABLE
 //fluid: not compressable (water in a fallon)
// compressible: air - can get more dense
// this code only works for incompressible fluids

// VECTOR FIELD
// fluid will live inside of a box, a square actually. a grid of 512 by 512
// each has a velocity vector. if velocity vector is 0 then it will be still
// velocity field has xs and ys. can add z for third dimension
//
// DYE
// DENSITY OF DYE, THING TO ADD TO VISUALIZE. IMAGING PUTTING DYE ONT HE GRID
// need an array to store all xs , ys, and amount of dye for each stop.
// need previous and next of each x,y, and dye. so need 3
//
// need x,y, density and x0, y0 and density0

//
console.log('debug index');
var N = 128;
var iter = 16;
var SCALE = 4;
var t = 0;
var fluid;

function setup() {
  createCanvas(N*SCALE,N*SCALE);
  fluid = new fluidCLass(0.2, 0, 0.0000001);

}

function draw() {
  background(0);
  var cx = (0.5*width)/SCALE;
  var cy = (0.5*height)/SCALE;
  for (var i = -1; i <= 1; i++) {
    for (var j = -1; j <= 1; j++) {
      fluid.addDensity(cx+i, cy+j, random(50, 150));
    }
  }
  for (var i = 0; i < 2; i++) {
    var angle = noise(t) * TWO_PI * 2;
    var v = p5.Vector.fromAngle(angle);
    v.mult(0.2);
    t += 0.01;
    fluid.addVelocity(cx, cy, v.x, v.y );
  }

  fluid.step();
  fluid.renderD();
}

// function mouseDragged() {
//   //add dye
//   fluid.addDensity(mouseX/SCALE, mouseY/SCALE, 100);
//   fluid.renderD();
// }
