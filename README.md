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
[ShinyProxy](www.shinyproxy.io) to allow multiple users to access workHORSE app
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

## Population data

workHORSE scales its synthetic population to ONS population counts for every
year from 2003 to 2047, so simulations can run up to 2047. Local authorities use
April 2023 boundaries (296 in England).

- Up to 2025: ONS mid-year estimates.
- 2026–2047, regions and local authorities: ONS 2022-based subnational
  population projections (migration category variant, which ONS uses in place
  of a principal projection for this release).
- 2026–2047, England: ONS 2024-based national population projections
  (principal projection). From 2026, England is therefore not the sum of its
  regions or local authorities.

Sources, licence and how to regenerate the files are in
[`ONS_data/pop_size/README.md`](ONS_data/pop_size/README.md).

## Health economics: discounting and prices

The Output tab discounts costs and QALYs before it calculates cumulative
values, net monetary benefit (NMB), incremental cost-effectiveness ratios
(ICERs), benefit:cost ratios and the inequality indices. The settings are in the
Health Economics menu of the Output tab, where each has an information icon. The
explanatory text on the Dashboard and Cost-effectiveness tabs states the years
covered, the price year and the discount rates.

- **Year 0 is the first simulated year.** It is read from the simulation
  results, not from the Period slider, so moving that slider after a run does
  not change how results are discounted. Nothing is discounted in year 0. This
  follows the HM Treasury [Green
  Book](https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government)
  (2026, paragraph 6.10 and Table 8) and Annex A of its [supplementary guidance
  on
  discounting](https://www.gov.uk/government/publications/green-book-supplementary-guidance-discounting).
  The NICE [health technology evaluations
  manual](https://www.nice.org.uk/process/pmg36) (PMG36, 4.5.1) requires present
  values over the time horizon of the analysis but does not define year 0.
- **Formula.** A cost or a QALY that occurs t years after year 0 is multiplied
  by `1 / (1 + r)^t`, where r is the annual rate. At 3.5% the factors for years
  1 to 3 are 0.9662, 0.9335 and 0.9019, and at 1.5% they are 0.9852, 0.9707 and
  0.9563 (Annex A, Tables A.1 and A.2). The factors are applied to the annual
  values, before the cumulative sums.
- **Defaults.** 3.5% a year for costs, 1.5% a year for QALYs, and willingness to
  pay of £20,000 per QALY. The Green Book discounts costs at 3.5% and health
  effects at 1.5% in years 1 to 30, and at lower rates after that (3.0% and
  1.286% in years 31 to 75, Table 3.A of the supplementary guidance). workHORSE
  applies the rate you choose to every year, which makes no difference within 30
  years of year 0. The NICE reference case discounts costs and health effects at
  3.5%, with 1.5% for both as an alternative analysis in specific circumstances
  (PMG36, 4.5.1 to 4.5.3); to use the NICE rates, set both sliders to 3.5%.
  NICE guidelines (PMG20) generally consider an ICER below £20,000 per QALY
  gained cost effective, and technology appraisals (PMG36) have used £25,000 to
  £35,000 per QALY gained since April 2026.
- **What is discounted.** Every cost column (names ending in `_cost`, including
  the net and total costs) at the cost rate, and `net_utility` and total QALYs
  (`eq5d`) at the QALY rate. The "Most effective" ranking therefore compares
  discounted QALYs, and the relative inequality index divides discounted net
  QALYs by discounted total QALYs.
- **Prices.** Costs are in 2019 prices, the price year of the unit costs in
  `simulation/health_econ/input/` (as in O'Flaherty et al., *Health Technology
  Assessment* 2021;25(35)), so enter scenario costs in 2019 prices. Costs are
  discounted but not uprated, whereas the Green Book (paragraph 6.47) would
  express them in year 0 prices.
- **CSV downloads.** The summarised results and the raw model output end with
  six columns that record the settings: `discount_year0`,
  `discount_rate_costs_pct`, `discount_rate_qalys_pct`, `price_year`,
  `perspective` and `wtp_gbp_per_qaly`.

### Change from the previous method

Costs and QALYs in year `y` used to be multiplied by `(1 - r)^(y - 2019)`, which
discounts to 2019, as in the HTA report cited above. Two things have changed:

1. The factor is now the standard `1 / (1 + r)^t`. The previous multiplier was
   smaller by `(1 - r^2)^t`, which is 1.2% after ten years at 3.5%, so it
   over-discounted.
2. Year 0 is now the first simulated year instead of 2019. Costs and QALYs are
   discounted at different rates, so with 2019 as year 0 their relative weights,
   and therefore NMB and ICERs, depended on how far the first simulated year was
   from 2019. For a cost in 2031 in a simulation that starts in 2021, the
   multiplier at 3.5% is now 0.7089 (it was 0.6521). For a QALY at 1.5% it is
   now 0.8617 (it was 0.8341).

For example, take a simulation of Rutland from 2021 to 2047 with ten Monte Carlo
iterations, comparing a scenario that adds a smoking cessation programme to
health checks against the baseline scenario, at the default settings and with
the healthcare perspective. With the new method the median cumulative net QALYs
are 162.4 instead of 156.9, the incremental cost-effectiveness ratio in the
explanatory text is £10,851 instead of £10,343 per QALY (4.9% higher), and the
benefit:cost ratio is 1.89 instead of 1.98. The size of the change depends on
the scenario and on how far the first simulated year is from 2019.

### Differential discounting

With different rates for costs and QALYs, an ICER depends on how long after year
0 a scenario's costs and QALYs occur (Keeler and Cretin, *Management Science*
1983;29(3):300-306). If the whole stream of a scenario starts d years later, its
discounted costs are divided by `1.035^d` and its discounted QALYs by `1.015^d`,
so its ICER is multiplied by `(1.015 / 1.035)^d` through discounting alone: by
0.981 for d = 1, 0.962 for d = 2 and 0.907 for d = 5. A scenario that starts
later therefore looks more cost effective only because of discounting. Take care
when comparing scenarios whose programmes start in different years, or use the
same rate for costs and QALYs.
