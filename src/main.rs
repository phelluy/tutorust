const L: f64 = 1.;

fn peak(x: f64) -> f64 {
    let r2 = x * x;
    let eps = 0.2;
    let eps2 = eps * eps;
    if r2 / eps2 < 1. {
        (1. - r2 / eps2).powi(4)
    } else {
        0.
    }
}

fn bosse2d(x: f64, y: f64) -> f64 {
    let r = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)).sqrt();
    peak(r)
}

fn bosse(x: f64) -> f64 {
    peak(x-0.5)
}

fn writepycom() {
    // use std::path::Path;
    // if 
    let pycom = r#"
from matplotlib.pyplot import *
from math import *
import numpy as np
with open("skyrs_plotpy.dat", "r") as f:
    contenu = f.read().split("\n\n")
    # print(contenu)
x = contenu[0].split()
nx = len(x) - 1
x = np.array([float(x[i]) for i in range(nx+1)])
# print(x)
y = contenu[1].split()
ny = len(y) -1
y = np.array([float(y[i]) for i in range(ny+1)])
# print(y)
x , y = np.meshgrid(x,y)
z = contenu[2].split()
nz = len(z)
z = np.array([float(z[i]) for i in range(nz)]).reshape((ny+1,nx+1))
# print(z)
# plot(xi, unum, color="blue")
# plot(xi, uex, color="red")
# def quit_figure(event):
#     if event.key == 'q':
#         close(event.canvas.figure)
# cid = gcf().canvas.mpl_connect('key_press_event', quit_figure)
fig, ax = subplots()
ax.axis('equal')
cs = ax.pcolormesh(x,y,z,cmap = 'jet')
# cs = ax.countourf(x,y,z,100,cmap = 'jet')
cbar = fig.colorbar(cs)
print("press \'q\' to quit...");
show()
"#;
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;

    let meshfile = File::create("skyrs_plot.py").unwrap();
    let mut meshfile = BufWriter::new(meshfile); // create a buffer for faster writes...
    writeln!(meshfile, "{}", pycom).unwrap();
}
/// Plots a 2D data set using matplotlib.
fn plotpy(xp: &Vec<f64>, yp: &Vec<f64>, zp: &Vec<f64>) {
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    {
        writepycom();
        let meshfile = File::create("skyrs_plotpy.dat").unwrap();
        let mut meshfile = BufWriter::new(meshfile); // create a buffer for faster writes...
        xp.iter().for_each(|x| {
            writeln!(meshfile, "{}", x).unwrap();
        });
        writeln!(meshfile).unwrap();
        yp.iter().for_each(|y| {
            writeln!(meshfile, "{}", y).unwrap();
        });
        writeln!(meshfile).unwrap();
        zp.iter().for_each(|z| {
            writeln!(meshfile, "{}", z).unwrap();
        });
    } // ensures that the file is closed

    use std::process::Command;
    Command::new("python3")
        .arg("skyrs_plot.py")
        .status()
        .expect("Plot failed: you need Python3 and Matplotlib in your PATH.");
}


#[derive(Debug, Clone)]
struct Data2D<T> {
    n: usize,
    data: Vec<T>,
}

use std::ops::Index;
impl<T> Index<usize> for Data2D<T> {
    type Output = [T];
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index * self.n .. (index+1) * self.n]
    }
}
// voir https://stackoverflow.com/questions/33770989/implementing-the-index-operator-for-matrices-with-multiple-parameters
use std::ops::IndexMut;
impl<T> IndexMut<usize> for Data2D<T> {
    //type Output = [T];
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index * self.n .. (index+1) * self.n]
    }
}

impl<T> Data2D<T> {
    fn new(n: usize, init_fn: fn(f64,f64)->T) -> Data2D<T> {
        let dx = L / n as f64;
        let data = (0..n*n).map(|k| {
            let i = k / n;
            let j = k%n;
            init_fn(i as f64 * dx, j as f64 * dx)
        }).collect();
        Data2D {
            n,
            data,
        }
    }
}
impl Data2D<f64> {
    fn plot(&self) {
        let dx = L/self.n as f64;
        let xc: Vec<f64> = (0..self.n).map(|i| (i as f64 -0.5) * dx ).collect();
        let yc = xc.clone();
        plotpy(&xc,&yc,&self.data)
    }
}


fn main() {
    let nx = 10;

    let mut tab2d = Data2D::new(nx,bosse2d);

    println!("data={}",tab2d[5][5]);

    tab2d[0][0] = 10.;

    tab2d.plot();

    let dx = L / nx as f64;

    let tmax = 3.;

    let cfl = 0.8;

    let dt = dx * cfl;

    let mut u_now = vec![0.; nx + 2];
    let mut unext = vec![0.; nx + 2];
    let mut uprev = vec![0.; nx + 2];

    let mut xc = vec![0.; nx + 2];

    for i in 0..nx + 2 {
        xc[i] = i as f64 * dx - dx / 2.;
    }

    for i in 0..nx + 2 {
        u_now[i] = bosse(xc[i]);
        uprev[i] = bosse(xc[i]);
    }

    let mut t = 0.;

    let b = dt / dx;
    let mut count = 0;
    let plotfreq = 100;

    while t < tmax {
        if count % plotfreq == 0 {
            plot1d(&xc, &u_now, &u_now);
        }
        count += 1;

        for i in 1..nx + 1 {
            unext[i] =
                -uprev[i] + 2. * (1. - b * b) * u_now[i] + b * b * (u_now[i - 1] + u_now[i + 1]);
        }
        //unext[0] = unext[1];
        // fixed the left 
        unext[0] = 0.;
        // free on the right
        unext[nx + 1] = unext[nx];
        //unext[nx + 1] = 0.;
        for i in 0..nx + 2 {
            uprev[i] = u_now[i];
            u_now[i] = unext[i];
        }
        t += dt;

        println!("t={}, dt={}", t, dt);
    }

    plot1d(&xc, &u_now, &u_now);
}

fn plot1d(x: &Vec<f64>, y: &Vec<f64>, z: &Vec<f64>) {
    let filename = "ploplo.dat";
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    {
        let meshfile = File::create(filename).unwrap();
        let mut meshfile = BufWriter::new(meshfile); // create a buffer for faster writes...

        x.iter()
            .zip(y.iter().zip(z.iter()))
            .for_each(|(x, (y, z))| {
                writeln!(meshfile, "{} {} {}", *x, *y, *z).unwrap();
            });
    }

    use std::process::Command;
    Command::new("python3")
        .arg("src/plot1d.py")
        .status()
        .expect("plot failed !");
}

