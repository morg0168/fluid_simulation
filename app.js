
console.log('debug app');
//lookup function for grid, returns one dimensional index -2 location from 1 d arrya
function IX(a,b) {
   var x = constrain(a, 0, N-1);
   var y = constrain(b, 0, N-1);
  return x + (y * N);
}

var fluidCLass = function(dt, diffusion, viscosity) {
  this.size = N;
  this.dt = dt; //time step
  this.diff = diffusion; //controls the velocity and devort dye diffuses throughout fluid
  this.visc = viscosity; //thickness of fluid

  this.s = new Array(N*N).fill(0);
  this.density = new Array(N*N).fill(0);

  this.Vx = new Array(N*N).fill(0);
  this.Vy = new Array(N*N).fill(0);

  this.Vx0 = new Array(N*N).fill(0);
  this.Vy0 =new Array(N*N).fill(0);

  this.step = function() {

    diffuse(1, this.Vx0, this.Vx, this.visc, this.dt);
    diffuse(2, this.Vy0, this.Vy, this.visc, this.dt);

    project(this.Vx0, this.Vy0, this.Vx, this.Vy);

    advect(1, this.Vx, this.Vx0, this.Vx0, this.Vy0, this.dt);
    advect(2, this.Vy, this.Vy0, this.Vx0, this.Vy0, this.dt);

    project(this.Vx, this.Vy, this.Vx0, this.Vy0);

    diffuse(0, this.s, this.density, this.diff, this.dt);
    advect(0, this.density, this.s, this.Vx, this.Vy, this.dt);
  };
  this.addDensity = function (x, y, amount) {
    var index = IX(x,y);
    this.density[index] += amount;
  };
  this.addVelocity= function(x, y, amountX, amountY) {
    var index = IX(x,y);
    this.Vx[index] += amountX;
    this.Vy[index] += amountY;
  };
  this.renderD = function() {
    colorMode(HSB, 255);
    for (var i = 0; i < N; i++){
      for (var j = 0; j < N; j++){
        var x = i * SCALE;
        var y = j * SCALE;
         var d = this.density[IX(i,j)];
         fill((d + 50) % 255, 200, d);
         noStroke();
         rect(x,y, SCALE, SCALE);
      }
    }
  };
  this.renderV = function() {
    for (var i = 0; i < N; i++){
      for (var j = 0; j < N; j++){
        var x = i * SCALE;
        var y = j * SCALE;
        var vx = this.Vx[IX(i,j)];
        var vy = this.Vy[IX(i,j)];
        stroke(255);

        if (!(Math.abs(vx) < 0.1 && Math.abs(vy) <= 0.1)) {
          line(x,y,x+vx*SCALE, y+vy*SCALE);
        }
      }
    }
  };
  this.fadeD = function() {
    for (var i =0; i < this.density.length; i++) {
      var d = density[i];
      this.density[i] = constrain(d-0.02, 0, 255);
    }
  };
}


//can diffuse any arbitary array of numbers versus previous array of numbers
//spreading out
// 1, this.Vx0, this.Vx, this.visc, this.dt
function diffuse(b,x, x0, diff, dt) {
  var a = dt * diff * (N-2) * (N-2);
  lin_solve(b,x,x0,a,1+6 * a);
}

function lin_solve(b,x,x0,a, c) {
  var cRecip = 1.0 /c;
  //New value of cell is a function of itself and all of its neighbours
  for (var k = 0; k < iter; k++) {
    for (var j = 1; j < N - 1; j++) {
      for (var i = 1; i < N -1; i++) {
       x[IX(i, j)] =
           (x0[IX(i, j)]
           + a*(    x[IX(i+1, j)]
           +x[IX(i-1, j)]
           +x[IX(i, j+1)]
           +x[IX(i, j-1)]
           )) * cRecip;
      }
    }
    set_bnd(b,x);
  }
}

function project (velocX, velocY, p, div) {
  var velocX = velocX;
  var velocY = velocY;
  var p = p;
  var div = div;
    for (var j = 1; j < N - 1; j++) {
      for (var i = 1; i < N -1; i++) {
         div[IX(i,j)] = -0.5*(
          + velocX[IX(i+1,j)]
          - velocX[IX(i-1,j)]
          + velocY[IX(i,j + 1)]
          - velocY[IX(i,j-1)]
        )/N;
         p[IX(i,j)] = 0;
      }
    }
    set_bnd(0, div);
    set_bnd(0,p);
    lin_solve(0,p,div,1,6);
    for (var j = 1; j < N - 1; j++) {
      for (var i = 1; i < N -1; i++) {
         velocX[IX(i, j)] -= 0.5 * (  p[IX(i+1, j)]
         -p[IX(i-1, j)]) * N;
         velocY[IX(i, j)] -= 0.5 * (  p[IX(i, j+1)]
         -p[IX(i, j-1)]) * N;
      }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

//advection - motion associated with the addVelocity
function advect(b, d, d0, velocX, velocY, dt) {
  var i0, i1, j0, j1;

  var dtx = dt * (N - 2);
  var dty = dt * (N - 2);

  var s0, s1, t0, t1;
  var tmp1, tmp2, x, y;

  var Nfloat = N;
  var ifloat, jfloat;
  var i, j;

  for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
    for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
      tmp1 = dtx * velocX[IX(i, j)];
      tmp2 = dty * velocY[IX(i, j)];
      x    = ifloat - tmp1;
      y    = jfloat - tmp2;

      if (x < 0.5) {
        x = 0.5;
      }
      if (x > Nfloat + 0.5) {
        x = Nfloat + 0.5;
      }
      i0 = Math.floor(x);
      i1 = i0 + 1.0;
      if (y < 0.5) {
        y = 0.5;
      }
      if (y > Nfloat + 0.5) {
        y = Nfloat + 0.5;
      }
      j0 = Math.floor(y);
      j1 = j0 + 1.0;

      s1 = x - i0;
      s0 = 1.0 - s1;
      t1 = y - j0;
      t0 = 1.0 - t1;

      var i0i = parseInt(i0);
      var i1i = parseInt(i1);
      var j0i = parseInt(j0);
      var j1i = parseInt(j1);

      // DOUBLE CHECK THIS!!!
      d[IX(i, j)] =
        (s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
        s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]));
    }
  }

  set_bnd(b, d);
}

//keep fluid in box
//b= which wall
//every velocity in the last reverses velocity in edge columns or rows according the next to column or rows
//corners does average of two neighbours
//wall gets velicty that perfectly counters it
function set_bnd (b, x) {
  for (var i = 1; i < N - 1; i++) {
    x[IX(i, 0  )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1 )];
    x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
  }
  for (var j = 1; j < N - 1; j++) {
    x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
    x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
  }

  x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
  x[IX(0, N-1)] = 0.5 * (x[IX(1, N-1)] + x[IX(0, N-2)]);
  x[IX(N-1, 0)] = 0.5 * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
  x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);
}
