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
    let r = ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)).sqrt();
    peak(r)
}

fn writepycom() {
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
#cs = ax.pcolormesh(x,y,z,cmap = 'jet',vmin=-0.2,vmax=0.2)
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
        &self.data[index * self.n..(index + 1) * self.n]
    }
}

use std::ops::IndexMut;
impl<T> IndexMut<usize> for Data2D<T> {
    //type Output = [T];
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index * self.n..(index + 1) * self.n]
    }
}

impl<T> Data2D<T> {
    fn new(n: usize, init_fn: fn(f64, f64) -> T) -> Data2D<T> {
        let dx = L / (n - 2) as f64;
        let data = (0..n * n)
            .map(|k| {
                let i = k / n;
                let j = k % n;
                init_fn(i as f64 * dx - dx / 2., j as f64 * dx - dx / 2.)
            })
            .collect();
        Data2D { n, data }
    }
}

impl Data2D<f64> {
    fn plot(&self) {
        let dx = L / (self.n - 2) as f64;
        let xc: Vec<f64> = (0..self.n).map(|i| (i as f64 - 0.5) * dx).collect();
        let yc = xc.clone();
        plotpy(&xc, &yc, &self.data)
    }
}

fn main() {
    let nx = 1000;

    let dx = L / nx as f64;

    let tmax = 0.6;

    let cfl = 0.4;

    let dt = dx * cfl;

    //let xc: Vec<f64> = (0..nx + 2).map(|i| i as f64 * dx - dx / 2.).collect();

    let mut u_now = Data2D::new(nx + 2, bosse2d);
    let mut unext = u_now.clone();
    let mut uprev = u_now.clone();

    let mut t = 0.;

    let b = dt / dx;
    let mut count = 0;
    let plotfreq = 200000;

    while t < tmax {
        if count % plotfreq == 0 {
            u_now.plot();
        }
        count += 1;
        for i in 1..nx + 1 {
            for j in 1..nx + 1 {
                unext[i][j] = -uprev[i][j]
                    + 2. * (1. - 2. * b * b) * u_now[i][j]
                    + b * b * (u_now[i - 1][j] + u_now[i + 1][j])
                    + b * b * (u_now[i][j - 1] + u_now[i][j + 1]);
            }
        }
        for k in 1..nx + 1 {
            unext[k][0] = unext[k][1];
            unext[k][nx + 1] = unext[k][nx];
            unext[0][k] = unext[1][k];
            unext[nx + 1][k] = unext[nx][k];
            // unext[k][0] = 0.;
            // unext[k][nx+1] = 0.;
            // unext[0][k] = 0.;
            // unext[nx+1][k] = 0.;
        }
        for i in 0..nx + 2 {
            for j in 0..nx + 2 {
                uprev[i][j] = u_now[i][j];
                u_now[i][j] = unext[i][j];
            }
        }
        t += dt;

        println!("t={}, dt={}", t, dt);
    }
    u_now.plot();
}
