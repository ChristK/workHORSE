# workHORSE microsimulation

--------------------------------------------------------------------------------

workHORSE is an implementation of the IMPACTncd framework, developed by Chris
Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been funded
by NIHR HTA Project: 16/165/01 - workHORSE: Health Outcomes Research Simulation
Environment. The views expressed are those of the authors and not necessarily
those of the NHS, the NIHR or the Department of Health. The main purpose of
workHORSE is for in-silico experimentation with different forms of Health
Checks, including the [NHS Health Check
Programme](https://www.healthcheck.nhs.uk/), in England.

Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos

workHORSE is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along
with this program; if not, see <http://www.gnu.org/licenses/> or write to the
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.

## workHORSE deployment instructions

The easiest way to deploy the workHORSE app is using a [Docker
container](https://www.docker.com/resources/what-container). workHORSE requires
a workstation with at least 20-cores and 256Gb RAM per concurrent user to run.
Many cloud-computing providers can fulfil these requirements nowadays relatively
cheap. Please note that these deployment instructions are not considering
security

We offer two ways to deploy workHORSE app. The first, will allow only one user
to access the app (single-user deployment). The second, leverages
[ShinyProxy][www.shinyproxy.io] to allow multiple users to access workHORSE app
concurrently, independently of each other (multi-user deployment).

We provide deployment instruction for two operating systems: Windows 10 (Pro,
Enterprise, or Edu versions only) and Ubuntu Linux. In reality workHORSE app
exploits the power and flexibility of [Docker
containers](https://www.docker.com/products/container-runtime) and can be
deployed in any operating system that is supported by [Docker](www.docker.com).

### Single-user installation

#### Linux – Ubuntu 20.04 LTS

1.  Open terminal and install Docker. Detailed instruction for Docker
    installation can be found
    [here](https://docs.docker.com/engine/install/ubuntu/).

2.  Get docker image with workHORSE app (this may take some time depending of
    your Internet connection).

``` bash
sudo docker pull chriskypri/workhorse-app
```

1.  Create docker volume for storing synthetic population data.

``` bash
sudo docker volume create workhorse-volume
```

1.  Run docker image:

``` bash
sudo docker run --mount source=workhorse-volume,target=/mnt/storage_fast/synthpop -p 9898:9898 -it chriskypri/workhorse-app
```

1.  Now you should be able to run the WorkHORSE app by opening web browser and
    open address **localhost:9898**

#### Windows 10 (not Home Edition)

1.  Download Docker Desktop for Windows
    [here](https://www.docker.com/get-started)
2.  Run Docker Desktop Installer  
    **Do not** check the option "Use Windows containers instead of Linux…" (see
    picture below)

![](www/images/608cfcc15c090dc41bebcf3c1458570a.png?raw=true)

1.  Restart Windows

2.  Configure Docker: Click Docker icon in the messaging area of Windows Desktop
    and go to 'Settings'

![](www/images/d841060d88640ee1d5b7571a625dc764.png?raw=true)

In Resources -\> Advanced select at least 4 CPUs and at least 8GB of memory.

Ideally you should select 20 CPUs and 256Gb of RAM if these are available in
your machine. Otherwise, the simulations may take several hours to complete, or
crash unexpectedly.

![](www/images/b24d31b4ba8461c7b6ca2a0b3c7dc3e6.png?raw=true)

Then click 'Apply & Restart'

1.  Open windows terminal (i.e. Windows PowerShell – press win key + R, then
    type **powershell**)

2.  Run commands:

``` bash
docker pull chriskypri/workhorse-app
docker volume create workhorse-volume
docker run --mount source=workhorse-volume,target=/mnt/storage_fast/synthpop -p 9898:9898 -it chriskypri/workhorse-app
```

1.  The window like below should appear. Allow docker to communicate via network
    interface

![](www/images/5a8401c5b8c394a55654afb0ae66fe5c.png?raw=true)

1.  Now you should be able to run the WorkHORSE app by opening web browser and
    open address: **localhost:9898**

### Multi-user installation (Linux – Ubuntu 20.04 LTS)

1.  Open terminal and install Docker. Detailed instruction for Docker
    installation can be found
    [here](https://docs.docker.com/engine/install/ubuntu/).

2.  Get docker image with workHORSE app (this may take some time depending of
    your Internet connection).

``` bash
sudo docker pull chriskypri/workhorse-app
```

1.  Get docker image with shinyproxy.

``` bash
sudo docker pull chriskypri/workhorse-shinyproxy
```

1.  Create docker volume for storing synthetic population data.

``` bash
sudo docker volume create workhorse-volume
```

1.  Run command

``` bash
sudo docker network create sp-workhorse-net
```

1.  Run shinyproxy image

``` bash
sudo docker run -d -v /var/run/docker.sock:/var/run/docker.sock --net sp-example-net -p 8080:8080 chriskypri/workhorse-shinyproxy
```

1.  Now you should be able to run the WorkHORSE app by opening web browser and
    open address **localhost:8080**

## Cloning this Repo

You can clone this repository, however, workHORSE uses some large files that
cannot be uploaded to GitHub repo. These files are uploaded to GitHub releases.
After you clone this GitHub repo, please source the included R script
`gh_deploy.R` to download the additional large files.
