/* IMPACTncd: A decision support tool for primary prevention of NCDs
 Copyright (C) 2015  Chris Kypridemos
 
 IMPACTncd is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/>
 or write to the Free Software Foundation, Inc., 51 Franklin Street,
 Fifth Floor, Boston, MA 02110-1301  USA.*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <string>
//#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
List meansd(NumericVector x) {
  x = na_omit(x);
  double meanm = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    meanm += x[i];
  }
  
  meanm /= n;
  
  double variance = 0;
  for(int i = 0; i < n; i++) 
  {
    variance += (x[i]-meanm)*(x[i]-meanm)/(n-1);
  }
  
  return Rcpp::List::create(meanm, sqrt(variance));
}

//quantile implementation for default R method (type 7)
// [[Rcpp::export]]
NumericVector fquantile(NumericVector x, NumericVector probs, bool na_rm = true) {
  if (all(is_na(x))) 
  {
    NumericVector out(probs.size(), NA_REAL);
    return(out);
  }
  if (na_rm) x = na_omit(x);
  const int n = x.size();
  NumericVector out(probs.size());
  IntegerVector ii(probs.size());
  NumericVector h(probs.size());
  NumericVector index = 1 + (n - 1) * probs;
  NumericVector lo = floor(index); //floor
  //ceiling
  NumericVector hi = ceiling(index);
  for(int i = 0; i < probs.size(); i++) 
  {//catch corner case when index element is int and ceiling = floor
    h[i] = index[i] - lo[i];
  }
  std::sort(x.begin(), x.end());
  out = x[as<IntegerVector>(lo) - 1];
  if (all(is_na(ii))) 
  {
    return(out);
  } else
  {
    x = x[as<IntegerVector>(hi) - 1];
    for(int i = 0; i < probs.size(); i++) 
    {
      out[i] = (1 - h[i]) * out[i] + h[i] * x[i];
    }
    return(out);
  }
} 


// [[Rcpp::export]]
NumericMatrix fquantile_byid(NumericVector x, NumericVector q, IntegerVector id) {
  // Need to be sorted by id 
  const int n = x.size();
  NumericMatrix z(unique(id).size(), q.size());
  int start = 0;
  int counter = 0;
  int end = 0;
  int counter_row = 0;
  NumericVector par_out(q.size());
  
  for (int i = 1; i < n; i++) { // start from 2nd element
    if (id[i] == id[i-1]) counter++;
    else
    {
      start = i - counter - 1;
      end = i - 1;
      counter = 0;
      z.row(counter_row) = fquantile(x[seq(start, end)], q);
      counter_row++;
    }
  }
  // take care the last group
  z.row(counter_row) =  fquantile(x[seq(n-counter-1, n-1)], q);
  return(z);
}

// The following code was adapted from QDiabetes-2013 (http://qdiabetes.org, http://svn.clinrisk.co.uk/opensource/qdiabetes)
// This is about 25% slower than when input is vectors instead of DF
// [[Rcpp::export]]
DataFrame QDriskDF(DataFrame df, int surv) {
  IntegerVector age = df["age"];
  IntegerVector sex = df["sex"];
  IntegerVector b_corticosteroids = df["b_corticosteroids"];
  //IntegerVector b_cvd = df["b_cvd"];
  IntegerVector b_treatedhyp = df["bpmed"];
  NumericVector bmi = df["bmival.cvdlag"];
  IntegerVector ethrisk = df["origin"];
  IntegerVector fh_diab = df["fh.diab"];
  IntegerVector smoke_cat = df["smoke_cat"];
  NumericVector town = df["townsend"];
  
  const int n = df.nrows();
  NumericVector out_nocvd(n);
  NumericVector out_cvd(n);
  double a = 0.0;
  double b = 0.0;
  double dage = 0.0;
  double age_1 = 0.0;
  double age_2 = 0.0;
  double dbmi = 0.0;
  double bmi_1 = 0.0;
  double bmi_2 = 0.0;
  double dtown = 0.0;
  
  /* The conditional arrays */
  const double survivor_M[16] = {
    0.0,
    0.998213708400726,
    0.996353209018707,
    0.994382798671722,
    0.992213606834412,
    0.989733397960663,
    0.987064540386200,
    0.984254062175751,
    0.981255292892456,
    0.977990627288818,
    0.974455237388611,
    0.970843732357025,
    0.967315018177032,
    0.963437378406525,
    0.959633111953735,
    0.955690681934357
  };
  const double Iethrisk_M[10] = {
    0.0,
    0.0,
    1.2366090720913343000000000,
    1.4716746107789032000000000,
    1.8073235649498174000000000,
    1.2056055595936399000000000,
    0.6032369975938766100000000,
    0.9095436207452737300000000,
    0.9137604632927512900000000,
    0.7123719045990779500000000
  };
  const double Ismoke_M[5] = {
    0.0,
    0.1618238582395977700000000,
    0.1902020385619117000000000,
    0.3210636179312467100000000,
    0.4140001301797494600000000
  }; 
  const double survivor_F[16] = {
    0.0,
    0.998714804649353,
    0.997435748577118,
    0.996052920818329,
    0.994562506675720,
    0.992949724197388,
    0.991141080856323,
    0.989293158054352,
    0.987293541431427,
    0.985133886337280,
    0.982810735702515,
    0.980465650558472,
    0.978020071983337,
    0.975493073463440,
    0.972945988178253,
    0.970350146293640
  };
  
  /* The conditional arrays */
  
  const double Iethrisk_F[10] = {
    0.0,
    0.0,
    1.2672136244963337000000000,
    1.4277605208830098000000000,
    1.8624060798103199000000000,
    1.2379988338989651000000000,
    0.4709034172907677900000000,
    0.3476400901703160500000000,
    1.1587283467731935000000000,
    0.7335499325010315100000000
  };
  const double Ismoke_F[5] = {
    0.0,
    0.1012537024947505100000000,
    0.1915520564380613400000000,
    0.3091894136143333900000000,
    0.4646730392693820800000000
  };
  
  for (int i = 0; i < n; i++)
  {
    if (sex[i] == 1) 
    {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage = (double) age[i];
      dage /= 10.0;
      age_1 = log(dage);
      age_2 = pow(dage, 3);
      dbmi = bmi[i];
      if (dbmi < 20.0) dbmi = 20.0;
      if (dbmi > 40.0) dbmi = 40.0;
      dbmi /=10.0;
      bmi_1 = pow(dbmi, 2);
      bmi_2 = pow(dbmi, 3);
      
      /* Centring the continuous variables */
      
      age_1 += -1.496771812438965;
      age_2 += -89.149559020996094;
      bmi_1 += -6.832604885101318;
      bmi_2 += -17.859918594360352;
      dtown  = town[i] + 0.132148191332817;
      
      /* Start of Sum */
      a = 0.0;
      
      /* The conditional sums */
      
      a += Iethrisk_M[ethrisk[i]];
      a += Ismoke_M[smoke_cat[i]];
      
      /* Sum from continuous values */
      
      a += age_1 * 4.4205598323371680000000000;
      a += age_2 * -0.0041132238299394193000000;
      a += bmi_1 * 1.1169895991721528000000000;
      a += bmi_2 * -0.1793529530251269100000000;
      a += dtown * 0.0291530815903822650000000;
      
      /* Sum from boolean values */
      
      a += b_corticosteroids[i] * 0.2059811979905692400000000;
      //a += b_cvd[i] * 0.3914728454990503100000000; // a for nocvd
      a += b_treatedhyp[i] * 0.5010787979849035100000000;
      a += fh_diab[i] * 0.8385800403428993500000000;
      
      /* Sum from interaction terms */
      
      a += age_1 * bmi_1 * 0.5051031253768063500000000;
      a += age_1 * bmi_2 * -0.1375233635462656000000000;
      a += age_1 * fh_diab[i] * -1.1463560542602569000000000;
      a += age_2 * bmi_1 * -0.0015800686452772700000000;
      a += age_2 * bmi_2 * 0.0003394090057824062300000;
      a += age_2 * fh_diab[i] * 0.0018524160353981260000000;
      
      b = a + 0.3914728454990503100000000; // b for cvd
      /* Calculate the score itself (0 to 1)*/
      out_nocvd[i] = 1.0 - pow(survivor_M[surv], exp(a));
      out_cvd[i]   = 1.0 - pow(survivor_M[surv], exp(b));
    }
    
    else if (sex[i] == 2) {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage = (double) age[i];
      dage /= 10.0;
      age_1 = pow(dage, 0.5);
      age_2 = pow(dage, 3);
      dbmi = bmi[i];
      if (dbmi < 20.0) dbmi = 20.0;
      if (dbmi > 40.0) dbmi = 40.0;
      dbmi /= 10.0;
      bmi_1 = dbmi;
      bmi_2 = pow(dbmi, 3);
      
      /* Centring the continuous variables */
      
      age_1 += -2.135220289230347;
      age_2 += -94.766799926757813;
      bmi_1 += -2.549620866775513;
      bmi_2 += -16.573980331420898;
      dtown  = town[i] + 0.224075347185135;
      
      /* Start of Sum */
      a = 0.0;
      
      /* The conditional sums */
      
      a += Iethrisk_F[ethrisk[i]];
      a += Ismoke_F[smoke_cat[i]];
      
      /* Sum from continuous values */
      
      a += age_1 * 4.3848331212989669000000000;
      a += age_2 * -0.0049763964406541149000000;
      a += bmi_1 * 3.3753336326064329000000000;
      a += bmi_2 * -0.0631628488667318330000000;
      a += dtown  * 0.0432726992998635970000000;
      
      /* Sum from boolean values */
      
      a += b_corticosteroids[i] * 0.2681990966241487000000000;
      //a += b_cvd * 0.3596176830984252900000000;
      a += b_treatedhyp[i] * 0.5314598436974725700000000;
      a += fh_diab[i] * 0.7315358845837640600000000;
      
      /* Sum from interaction terms */
      
      a += age_1 * bmi_1 * 1.3037832873997990000000000;
      a += age_1 * bmi_2 * -0.0708293717769046120000000;
      a += age_1 * fh_diab[i] * -0.7968266815834251800000000;
      a += age_2 * bmi_1 * -0.0067725323761278549000000;
      a += age_2 * bmi_2 * 0.0002374980728666116700000;
      a += age_2 * fh_diab[i] * 0.0017048228889394394000000;
      
      b = a + 0.3596176830984252900000000;
      /* Calculate the score itself */
      out_nocvd[i] = 1.0 - pow(survivor_F[surv], exp(a));
      out_cvd[i] = 1.0 - pow(survivor_F[surv], exp(b));
    }
  }
  return DataFrame::create(_["QDrisk_nocvd"]= out_nocvd, 
                           _["QDrisk_cvd"]= out_cvd);
}

// The following code was adapted from QDiabetes-2013 (http://qdiabetes.org, http://svn.clinrisk.co.uk/opensource/qdiabetes)
// [[Rcpp::export]]
List QDrisk(IntegerVector age, IntegerVector sex, IntegerVector b_corticosteroids, IntegerVector b_treatedhyp, 
             NumericVector bmi, IntegerVector ethrisk, IntegerVector fh_diab, IntegerVector smoke_cat, NumericVector town, 
             int surv) {
  
  
  const int n = age.size();
  NumericVector out_nocvd(n);
  NumericVector out_cvd(n);
  double a = 0.0;
  double b = 0.0;
  double dage = 0.0;
  double age_1 = 0.0;
  double age_2 = 0.0;
  double dbmi = 0.0;
  double bmi_1 = 0.0;
  double bmi_2 = 0.0;
  double dtown = 0.0;
  
  /* The conditional arrays */
  const double survivor_M[16] = {
    0.0,
    0.998213708400726,
    0.996353209018707,
    0.994382798671722,
    0.992213606834412,
    0.989733397960663,
    0.987064540386200,
    0.984254062175751,
    0.981255292892456,
    0.977990627288818,
    0.974455237388611,
    0.970843732357025,
    0.967315018177032,
    0.963437378406525,
    0.959633111953735,
    0.955690681934357
  };
  const double Iethrisk_M[10] = {
    0.0,
    0.0,
    1.2366090720913343000000000,
    1.4716746107789032000000000,
    1.8073235649498174000000000,
    1.2056055595936399000000000,
    0.6032369975938766100000000,
    0.9095436207452737300000000,
    0.9137604632927512900000000,
    0.7123719045990779500000000
  };
  const double Ismoke_M[5] = {
    0.0,
    0.1618238582395977700000000,
    0.1902020385619117000000000,
    0.3210636179312467100000000,
    0.4140001301797494600000000
  }; 
  const double survivor_F[16] = {
    0.0,
    0.998714804649353,
    0.997435748577118,
    0.996052920818329,
    0.994562506675720,
    0.992949724197388,
    0.991141080856323,
    0.989293158054352,
    0.987293541431427,
    0.985133886337280,
    0.982810735702515,
    0.980465650558472,
    0.978020071983337,
    0.975493073463440,
    0.972945988178253,
    0.970350146293640
  };
  
  /* The conditional arrays */
  
  const double Iethrisk_F[10] = {
    0.0,
    0.0,
    1.2672136244963337000000000,
    1.4277605208830098000000000,
    1.8624060798103199000000000,
    1.2379988338989651000000000,
    0.4709034172907677900000000,
    0.3476400901703160500000000,
    1.1587283467731935000000000,
    0.7335499325010315100000000
  };
  const double Ismoke_F[5] = {
    0.0,
    0.1012537024947505100000000,
    0.1915520564380613400000000,
    0.3091894136143333900000000,
    0.4646730392693820800000000
  };
  
  for (int i = 0; i < n; i++)
  {
    if (sex[i] == 1) 
    {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage = (double) age[i];
      dage /= 10.0;
      age_1 = log(dage);
      age_2 = pow(dage, 3);
      dbmi = bmi[i];
      if (dbmi < 20.0) dbmi = 20.0;
      if (dbmi > 40.0) dbmi = 40.0;
      dbmi /=10.0;
      bmi_1 = pow(dbmi, 2);
      bmi_2 = pow(dbmi, 3);
      
      /* Centring the continuous variables */
      
      age_1 += -1.496771812438965;
      age_2 += -89.149559020996094;
      bmi_1 += -6.832604885101318;
      bmi_2 += -17.859918594360352;
      dtown  = town[i] + 0.132148191332817;
      
      /* Start of Sum */
      a = 0.0;
      
      /* The conditional sums */
      
      a += Iethrisk_M[ethrisk[i]];
      a += Ismoke_M[smoke_cat[i]];
      
      /* Sum from continuous values */
      
      a += age_1 * 4.4205598323371680000000000;
      a += age_2 * -0.0041132238299394193000000;
      a += bmi_1 * 1.1169895991721528000000000;
      a += bmi_2 * -0.1793529530251269100000000;
      a += dtown * 0.0291530815903822650000000;
      
      /* Sum from boolean values */
      
      a += b_corticosteroids[i] * 0.2059811979905692400000000;
      //a += b_cvd[i] * 0.3914728454990503100000000; // a for nocvd
      a += b_treatedhyp[i] * 0.5010787979849035100000000;
      a += fh_diab[i] * 0.8385800403428993500000000;
      
      /* Sum from interaction terms */
      
      a += age_1 * bmi_1 * 0.5051031253768063500000000;
      a += age_1 * bmi_2 * -0.1375233635462656000000000;
      a += age_1 * fh_diab[i] * -1.1463560542602569000000000;
      a += age_2 * bmi_1 * -0.0015800686452772700000000;
      a += age_2 * bmi_2 * 0.0003394090057824062300000;
      a += age_2 * fh_diab[i] * 0.0018524160353981260000000;
      
      b = a + 0.3914728454990503100000000; // b for cvd
      /* Calculate the score itself (0 to 1)*/
      out_nocvd[i] = 1.0 - pow(survivor_M[surv], exp(a));
      out_cvd[i]   = 1.0 - pow(survivor_M[surv], exp(b));
    }
    
    else if (sex[i] == 2) {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage = (double) age[i];
      dage /= 10.0;
      age_1 = pow(dage, 0.5);
      age_2 = pow(dage, 3);
      dbmi = bmi[i];
      if (dbmi < 20.0) dbmi = 20.0;
      if (dbmi > 40.0) dbmi = 40.0;
      dbmi /= 10.0;
      bmi_1 = dbmi;
      bmi_2 = pow(dbmi, 3);
      
      /* Centring the continuous variables */
      
      age_1 += -2.135220289230347;
      age_2 += -94.766799926757813;
      bmi_1 += -2.549620866775513;
      bmi_2 += -16.573980331420898;
      dtown  = town[i] + 0.224075347185135;
      
      /* Start of Sum */
      a = 0.0;
      
      /* The conditional sums */
      
      a += Iethrisk_F[ethrisk[i]];
      a += Ismoke_F[smoke_cat[i]];
      
      /* Sum from continuous values */
      
      a += age_1 * 4.3848331212989669000000000;
      a += age_2 * -0.0049763964406541149000000;
      a += bmi_1 * 3.3753336326064329000000000;
      a += bmi_2 * -0.0631628488667318330000000;
      a += dtown  * 0.0432726992998635970000000;
      
      /* Sum from boolean values */
      
      a += b_corticosteroids[i] * 0.2681990966241487000000000;
      //a += b_cvd * 0.3596176830984252900000000;
      a += b_treatedhyp[i] * 0.5314598436974725700000000;
      a += fh_diab[i] * 0.7315358845837640600000000;
      
      /* Sum from interaction terms */
      
      a += age_1 * bmi_1 * 1.3037832873997990000000000;
      a += age_1 * bmi_2 * -0.0708293717769046120000000;
      a += age_1 * fh_diab[i] * -0.7968266815834251800000000;
      a += age_2 * bmi_1 * -0.0067725323761278549000000;
      a += age_2 * bmi_2 * 0.0002374980728666116700000;
      a += age_2 * fh_diab[i] * 0.0017048228889394394000000;
      
      b = a + 0.3596176830984252900000000;
      /* Calculate the score itself */
      out_nocvd[i] = 1.0 - pow(survivor_F[surv], exp(a));
      out_cvd[i] = 1.0 - pow(survivor_F[surv], exp(b));
    }
  }
  return Rcpp::List::create(_["QDrisk_nocvd"]= out_nocvd,
                            _["QDrisk_cvd"]= out_cvd);
}

// The following code was adapted from QRISK2-2014 (http://qrisk.org, original sources at http://svn.clinrisk.co.uk/opensource/qrisk2).
// [[Rcpp::export]]
double cvd_raw(
    int age,             // 25:84
    int sex,             // 1, 2
    int b_AF,            // 0, 1
    int b_ra,            // 0, 1
    int b_renal,         // 0, 1
    int b_treatedhyp,    // 0, 1
    int b_type1,         // 0, 1
    int b_type2,         // 0, 1
    double bmi,          // 20:40
    int ethrisk,         // 1:9
    int fh_cvd,          // 0, 1
    double rati,         // 1 - 12
    double sbp,          // 70 - 210
    int smoke_cat,       // 0:4
    double town          // -7 to 11
)
{
  if (bmi < 20.0) bmi = 20.0;
  if (bmi > 40.0) bmi = 40.0;
  
  if (rati < 1.0) rati = 1.0;
  if (rati > 12.0) rati = 12.0;
  
  if (sbp < 70.0) sbp = 70.0;
  if (sbp > 210.0) sbp = 210.0;
  
  b_type2 -= 1; // diabtotr in my model is factor 1,2 and converted to int 1, 2. I need 0,1 
  
  if (sex == 1) {
    double survivor[16] = {
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.977699398994446,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    };
    
    /* The conditional arrays */
    
    double Iethrisk[10] = {
      0.0,
      0.0,
      0.3567133647493443400000000,
      0.5369559608176189800000000,
      0.5190878419529624300000000,
      0.2182992106490147000000000,
      -0.3474174705898491800000000,
      -0.3674730037922803700000000,
      -0.3749664891426142700000000,
      -0.1926947742531604500000000
    };
    double Ismoke[5] = {
      0.0,
      0.2784649664157046200000000,
      0.6067834395168959500000000,
      0.7103835060989258700000000,
      0.8626172339181202900000000
    };
    
    /* Applying the fractional polynomial transforms */
    /* (which includes scaling)                      */
    
    double dage = age;
    dage=dage/10.0;
    double age_1 = pow(dage,-1);
    double age_2 = pow(dage,2);
    double dbmi = bmi;
    dbmi=dbmi/10.0;
    double bmi_1 = pow(dbmi,-2);
    double bmi_2 = pow(dbmi,-2)*log(dbmi);
    
    /* Centring the continuous variables */
    
    age_1 = age_1 - 0.232008963823318;
    age_2 = age_2 - 18.577636718750000;
    bmi_1 = bmi_1 - 0.146408438682556;
    bmi_2 = bmi_2 - 0.140651300549507;
    rati  = rati - 4.377167701721191;
    sbp   = sbp - 131.038314819335940;
    town  = town - 0.151332527399063;
    
    /* Start of Sum */
    double a = 0.0;
    
    /* The conditional sums */
    
    a += Iethrisk[ethrisk];
    a += Ismoke[smoke_cat];
    
    /* Sum from continuous values */
    
    a += age_1 * -17.6225543381945610000000000;
    a += age_2 * 0.0241873189298273640000000;
    a += bmi_1 * 1.7320282704272665000000000;
    a += bmi_2 * -7.2311754066699754000000000;
    a += rati * 0.1751387974012235100000000;
    a += sbp * 0.0101676305179196900000000;
    a += town * 0.0298177271496720960000000;
    
    /* Sum from boolean values */
    
    a += b_AF * 0.9890997526189402300000000;
    a += b_ra * 0.2541886209118611200000000;
    a += b_renal * 0.7949789230438320000000000;
    a += b_treatedhyp * 0.6229359479868044100000000;
    a += b_type1 * 1.3330353321463930000000000;
    a += b_type2 * 0.9372956828151940400000000;
    a += fh_cvd * 0.5923353736582422900000000;
    
    /* Sum from interaction terms */
    
    a += age_1 * (smoke_cat==1) * 0.9243747443632776000000000;
    a += age_1 * (smoke_cat==2) * 1.9597527500081284000000000;
    a += age_1 * (smoke_cat==3) * 2.9993544847631153000000000;
    a += age_1 * (smoke_cat==4) * 5.0370735254768100000000000;
    a += age_1 * b_AF * 8.2354205455482727000000000;
    a += age_1 * b_renal * -3.9747389951976779000000000;
    a += age_1 * b_treatedhyp * 7.8737743159167728000000000;
    a += age_1 * b_type1 * 5.4238504414460937000000000;
    a += age_1 * b_type2 * 5.0624161806530141000000000;
    a += age_1 * bmi_1 * 33.5437525167394240000000000;
    a += age_1 * bmi_2 * -129.9766738257203800000000000;
    a += age_1 * fh_cvd * 1.9279963874659789000000000;
    a += age_1 * sbp * 0.0523440892175620200000000;
    a += age_1 * town * -0.1730588074963540200000000;
    a += age_2 * (smoke_cat==1) * -0.0034466074038854394000000;
    a += age_2 * (smoke_cat==2) * -0.0050703431499952954000000;
    a += age_2 * (smoke_cat==3) * 0.0003216059799916440800000;
    a += age_2 * (smoke_cat==4) * 0.0031312537144240087000000;
    a += age_2 * b_AF * 0.0073291937255039966000000;
    a += age_2 * b_renal * -0.0261557073286531780000000;
    a += age_2 * b_treatedhyp * 0.0085556382622618121000000;
    a += age_2 * b_type1 * 0.0020586479482670723000000;
    a += age_2 * b_type2 * -0.0002328590770854172900000;
    a += age_2 * bmi_1 * 0.0811847212080794990000000;
    a += age_2 * bmi_2 * -0.2558919068850948300000000;
    a += age_2 * fh_cvd * -0.0056729073729663406000000;
    a += age_2 * sbp * -0.0000536584257307299330000;
    a += age_2 * town * -0.0010763305052605857000000;
    
    /* Calculate the score itself */
    double score = 1.0 - pow(survivor[10], exp(a));
    return score;
  }
  
  else if (sex == 2) {
    double survivor[16] = {
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.988948762416840,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    };
    
    /* The conditional arrays */
    
    double Iethrisk[10] = {
      0.0,
      0.0,
      0.2671958047902151500000000,
      0.7147534261793343500000000,
      0.3702894474455115700000000,
      0.2073797362620235500000000,
      -0.1744149722741736900000000,
      -0.3271878654368842200000000,
      -0.2200617876129250500000000,
      -0.2090388032466696800000000
    };
    double Ismoke[5] = {
      0.0,
      0.1947480856528854800000000,
      0.6229400520450627500000000,
      0.7405819891143352600000000,
      0.9134392684576959600000000
    };
    
    /* Applying the fractional polynomial transforms */
    /* (which includes scaling)                      */
    
    double dage = age;
    dage=dage/10.0;
    double age_1 = pow(dage,.5);
    double age_2 = dage;
    double dbmi = bmi;
    dbmi=dbmi/10.0;
    double bmi_1 = pow(dbmi,-2);
    double bmi_2 = pow(dbmi,-2)*log(dbmi);
    
    /* Centring the continuous variables */
    
    age_1 = age_1 - 2.099778413772583;
    age_2 = age_2 - 4.409069538116455;
    bmi_1 = bmi_1 - 0.154046609997749;
    bmi_2 = bmi_2 - 0.144072100520134;
    rati  = rati - 3.554229259490967;
    sbp   = sbp - 125.773628234863280;
    town  = town - 0.032508373260498;
    
    /* Start of Sum */
    double a=0.0;
    
    /* The conditional sums */
    
    a += Iethrisk[ethrisk];
    a += Ismoke[smoke_cat];
    
    /* Sum from continuous values */
    
    a += age_1 * 3.8734583855051343000000000;
    a += age_2 * 0.1346634304478384600000000;
    a += bmi_1 * -0.1557872403333062600000000;
    a += bmi_2 * -3.7727795566691125000000000;
    a += rati * 0.1525695208919679600000000;
    a += sbp * 0.0132165300119653560000000;
    a += town * 0.0643647529864017080000000;
    
    /* Sum from boolean values */
    
    a += b_AF * 1.4235421148946676000000000;
    a += b_ra * 0.3021462511553648100000000;
    a += b_renal * 0.8614743039721416400000000;
    a += b_treatedhyp * 0.5889355458733703800000000;
    a += b_type1 * 1.6684783657502795000000000;
    a += b_type2 * 1.1350165062510138000000000;
    a += fh_cvd * 0.5133972775738673300000000;
    
    /* Sum from interaction terms */
    
    a += age_1 * (smoke_cat==1) * 0.6891139747579299000000000;
    a += age_1 * (smoke_cat==2) * 0.6942632802121626600000000;
    a += age_1 * (smoke_cat==3) * -1.6952388644218186000000000;
    a += age_1 * (smoke_cat==4) * -1.2150150940219255000000000;
    a += age_1 * b_AF * -3.5855215448190969000000000;
    a += age_1 * b_renal * -3.0766647922469192000000000;
    a += age_1 * b_treatedhyp * -4.0295302811880314000000000;
    a += age_1 * b_type1 * -0.3344110567405778600000000;
    a += age_1 * b_type2 * -3.3144806806620530000000000;
    a += age_1 * bmi_1 * -5.5933905797230006000000000;
    a += age_1 * bmi_2 * 64.3635572837688980000000000;
    a += age_1 * fh_cvd * 0.8605433761217157200000000;
    a += age_1 * sbp * -0.0509321154551188590000000;
    a += age_1 * town * 0.1518664540724453700000000;
    a += age_2 * (smoke_cat==1) * -0.1765395485882681500000000;
    a += age_2 * (smoke_cat==2) * -0.2323836483278573000000000;
    a += age_2 * (smoke_cat==3) * 0.2734395770551826300000000;
    a += age_2 * (smoke_cat==4) * 0.1432552287454152700000000;
    a += age_2 * b_AF * 0.4986871390807032200000000;
    a += age_2 * b_renal * 0.4393033615664938600000000;
    a += age_2 * b_treatedhyp * 0.6904385790303250200000000;
    a += age_2 * b_type1 * -0.1734316566060327700000000;
    a += age_2 * b_type2 * 0.4864930655867949500000000;
    a += age_2 * bmi_1 * 1.5223341309207974000000000;
    a += age_2 * bmi_2 * -12.7413436207964070000000000;
    a += age_2 * fh_cvd * -0.2756708481415109900000000;
    a += age_2 * sbp * 0.0073790750039744186000000;
    a += age_2 * town * -0.0487465462679640900000000;
    
    /* Calculate the score itself */
    double score = 1.0 - pow(survivor[10], exp(a)) ;
    return score;
  }
  else return 0;
}

// Vectorise above
// [[Rcpp::export]]
NumericVector QRisk(
    IntegerVector age,
    IntegerVector sex,
    IntegerVector b_AF,
    IntegerVector b_ra,
    IntegerVector b_renal,
    IntegerVector b_treatedhyp,
    IntegerVector b_type1,
    IntegerVector b_type2,
    NumericVector bmi,
    IntegerVector ethrisk,
    IntegerVector fh_cvd,
    NumericVector rati,
    NumericVector sbp,
    IntegerVector smoke_cat,
    NumericVector town
)
{
  int n = age.size();
  NumericVector out(n);
  
  for (int i = 0; i < n; i++) {
    out[i] = cvd_raw(age[i], sex[i], b_AF[i], b_ra[i],
                     b_renal[i], b_treatedhyp[i], b_type1[i],
                                                         b_type2[i], bmi[i], ethrisk[i], fh_cvd[i],
                                                                                               rati[i], sbp[i], smoke_cat[i], town[i]);
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector fbound(NumericVector x, float a, float b) {
  if (a > b) {float c = a; a = b; b = c;}; // ensure a < b
  const int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++) {
    if (Rcpp::NumericVector::is_na(x[i])) out[i] = NA_REAL;
    else
    {
      if (x[i] < a) out[i] = a;
      else if (x[i] > b) out[i] = b;
      else out[i] = x[i];
    } 
  }
  return out;
}

// [[Rcpp::export]]
NumericVector fbound_inplace(NumericVector x, float a, float b) {
  if (a > b) {float c = a; a = b; b = c;}; // ensure a < b
  const int n = x.size();
  for(int i = 0; i < n; i++) {
    if (x[i] < a) x[i] = a;
    else if (x[i] > b) x[i] = b;
  } 
  
  return x; // transforms input vector inplace as long as an numeric is passed to it. If an integer vector is passed to it from R, an impricit copy is happening so it is not inplace anymore
}

// [[Rcpp::export]]
IntegerVector fbound_inplace_int(IntegerVector x, int a, int b) {
  if (a > b) {float c = a; a = b; b = c;}; // ensure a < b
  const int n = x.size();
  for(int i = 0; i < n; i++) {
    if (x[i] < a) x[i] = a;
    else if (x[i] > b) x[i] = b;
  } 
  
  return x; // transforms input vector inplace as long as an numeric is passed to it. If an integer vector is passed to it from R, an impricit copy is happening so it is not inplace anymore
} 

// [[Rcpp::export]]
LogicalVector fequal(NumericVector x, double tol) {
  NumericVector y = na_omit(x);
  const int n = y.size();
  for (int i = 0; i < n; ++i) {
    if (y[i] - y[0] > tol || y[0] - y[i] > tol)
      return Rcpp::wrap(false);
  }
  
  return Rcpp::wrap(true);
}

// [[Rcpp::export]]
NumericVector fnormalise(NumericVector x) { // between 0, 1
  const int n = x.size();
  const double minx = min(na_omit(x));
  const double maxx = max(na_omit(x));
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = (x[i] - minx) / (maxx - minx);
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector fcompress(NumericVector x, float limit) {
  const int n = x.size();
  NumericVector out = clone(x);
  float sum_out = sum(na_omit(out));
  int counter = 0;
  while(pow(sum_out - limit, 2) > 0.1)
  {
    if (counter == 1000) break;
    sum_out = sum(na_omit(out));
    for(int i = 0; i < n; i++)
    {
      out[i] *= limit/sum_out;
      if (out[i] < 0.0) out[i] = 0.0;
      else if (out[i] > 1.0) out[i] = 1.0;
    }
    counter++;
  } 
  return out; // transforms input vector inplace
}

// [[Rcpp::export]]
int count_if(LogicalVector x) {
  const int n = x.size();
  int counter = 0;
  for(int i = 0; i < n; i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter;
}

// [[Rcpp::export]]
double prop_if(LogicalVector x, bool na_rm = true) {

  if (na_rm == true) x = na_omit(x); // remove NA from denominator
  const int n = x.size();
  int counter = 0;
  for(int i = 0; i < n; i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter/(double)n;
}

// [[Rcpp::export]]
NumericMatrix df2mat_numeric(DataFrame x) {
  int xsize=x.size();
  NumericMatrix y(x.nrows(),xsize);
  for(int i=0; i<xsize; i++){
    y(_,i) = NumericVector(x[i]);
  }
  y.attr("dimnames") = List::create( R_NilValue, x.attr("names") ) ;
  return(y);
}

// [[Rcpp::export]]
IntegerMatrix df2mat_integer(DataFrame x) {
  int xsize=x.size();
  IntegerMatrix y(x.nrows(),xsize);
  for(int i=0; i<xsize; i++){
    y(_,i) = IntegerVector(x[i]);
  }
  y.attr("dimnames") = List::create( R_NilValue, x.attr("names") ) ;
  return(y);
}


#include "utils.h"


template <typename T, typename U>
IntegerVector do_counts(const T&);

template <>
inline IntegerVector do_counts<IntegerVector, int>(const IntegerVector& x) {
  std::map< int, int, NACompare<int> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  return wrap(counts);
}

template <>
inline IntegerVector do_counts<LogicalVector, int>(const LogicalVector& x) {
  std::map< int, int, NACompare<int> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  IntegerVector output = wrap(counts);
  
  // yuck
  SEXP namesptr = Rf_getAttrib(output, R_NamesSymbol);
  for (int i=0; i < output.size(); ++i) {
    if (strcmp(CHAR(STRING_ELT(namesptr, i)), "0") == 0) {
      SET_STRING_ELT(namesptr, i, Rf_mkChar("FALSE"));
    }
    if (strcmp(CHAR(STRING_ELT(namesptr, i)), "1") == 0) {
      SET_STRING_ELT(namesptr, i, Rf_mkChar("TRUE"));
    }
  }
  
  return output;
  
}

template <>
inline IntegerVector do_counts<CharacterVector, SEXP>(const CharacterVector& x) {
  IntegerVector tmp = table(x);
  std::map< SEXP, int, NACompare<SEXP> > counts;
  for (int i=0; i < tmp.size(); ++i) {
    counts[ STRING_ELT( tmp.attr("names"), i ) ] = tmp[i];
  }
  IntegerVector output = wrap(counts);
  
  CharacterVector names = output.attr("names");
  CharacterVector::iterator it = std::find( names.begin(), names.end(), "NA" );
  if (it != names.end()) {
    *it = NA_STRING;
  }
  return output;
  
}

template <>
inline IntegerVector do_counts<NumericVector, double>(const NumericVector& x) {
  std::map< double, int, NACompare<double> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  IntegerVector output = wrap(counts);
  
  // explicitly use R's double-to-character coercion to get good names
  int m = counts.size();
  NumericVector keys = no_init(m);
  typedef std::map< double, int, NACompare<double> >::iterator MapItr;
  int i = 0;
  for (MapItr it = counts.begin(); it != counts.end(); ++it) {
    keys[i] = it->first;
    ++i;
  }
  CharacterVector names = Rf_coerceVector(keys, STRSXP);
  output.attr("names") = names;
  // fix names
  for (int i=0; i < output.size(); ++i) {
    bool samestr = strcmp(
      CHAR(STRING_ELT(output.attr("names"), i)),
      "-0"
    ) == 0;
    if (samestr) {
      SET_STRING_ELT(output.attr("names"), i, Rf_mkChar("0"));
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector tableRcpp(SEXP x) {
  switch (TYPEOF(x)) {
  case INTSXP: return table(as<IntegerVector>(x));
  case REALSXP: return table(as<NumericVector>(x));
  case STRSXP: return table(as<CharacterVector>(x));
  case LGLSXP: return table(as<LogicalVector>(x));
  default: {
    stop("unrecognized SEXP type");
    return R_NilValue;
  }
  }
}

// [[Rcpp::export]]
IntegerVector counts(SEXP x) {
  switch (TYPEOF(x)) {
  case REALSXP: return do_counts<NumericVector, double>(x);
  case STRSXP: return do_counts<CharacterVector, SEXP>(x);
  case INTSXP: return do_counts<IntegerVector, int>(x);
  case LGLSXP: return do_counts<LogicalVector, int>(x);
  default: {
    stop("unrecognized SEXP type");
    return R_NilValue;
  }
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame counts_DF(SEXP x) {
  Rcpp::DataFrame out = counts(x);
  return(out);
}

// [[Rcpp::export]]
NumericVector roll_mean(NumericVector dat,  int window) {
  const int n = dat.size();
  window =  (window > n ) ? n : window;
  
  NumericVector out(n);
  double summed = 0.0;
  for (int i=0; i < window; ++i) {
    summed += dat[i];
  }
  out[0] = summed / window;
  for (int i=window; i < n; ++i) {
    summed += dat[i] - dat[i-window];
    out[i-window+1] = summed / window;
    //fill
    if (i > (n-window) ) out[i] = dat[i];
  }
  return out;
}

// Dataframe input. Slower than separated vectors input
// [[Rcpp::export]]
DataFrame simsmokDF(DataFrame df, NumericMatrix pr_relap, int relap_cutoff) {
  
  //access the columns
  IntegerVector cigst1   = df["cigst1"]; 
  NumericVector pr_init  = df["pr_init"];
  NumericVector pr_cess  = df["pr_cess"];
  NumericVector dice     = df["dice"];
  IntegerVector id       = df["id"];
  IntegerVector age      = df["age"];
  IntegerVector sex      = df["sex"];
  IntegerVector qimd     = df["qimd"];
  IntegerVector endsmoke = df["endsmoke"];
  IntegerVector smokyrs  = df["smokyrs"];
  NumericVector cigdyal  = df["cigdyal"];
  NumericVector numsmok  = df["numsmok"];
  IntegerVector year     = df["year"];
  NumericVector cigdyal_rank = df["cigdyal.rank"];
  // id should be sorted by year
  const int n = df.nrows();
  int nrow = 0;
  
  if ((cigst1[0] == 1) && (dice[0] < pr_init[0])) { //non smoker
    cigst1[0] = 4;
    endsmoke[0] = 0;
    smokyrs[0] = smokyrs[0] + 1;
    numsmok[0] = 0;
    cigdyal[0] = 99; //to be calculated later
  }
  
  if (cigst1[0] == 4) { //current smoker
    if (dice[0] < pr_cess[0]) {
      cigst1[0] = 3;
      endsmoke[0] = 1;
      numsmok[0] = cigdyal[0];
      cigdyal[0] = 0;
    } else smokyrs[0] = smokyrs[0] + 1; 
  }
  
  if ((cigst1[0] == 3) && (endsmoke[0] < relap_cutoff)){
    switch (sex[0]) {
    case 1: nrow = qimd[0] - 1; break;
    case 2: nrow = qimd[0] + 4; break;
    }
    if (dice[0] < (pr_relap(nrow, endsmoke[0] + 1)))  {
      cigst1[0] = 4;
      endsmoke[0] = 0;
      smokyrs[0] = smokyrs[0] + 1;
      cigdyal[0] = numsmok[0];
      numsmok[0] = 0;
    } 
  }
  if (cigst1[0] == 3) endsmoke[0] = endsmoke[0] + 1;
  
  for (int i = 1; i < n; i++) // start from second element
  {
    if (id[i] == id[i-1]) { // if the same person, take account previous year
      if (cigst1[i-1] == 1) { // never smoker
        if (dice[i] < pr_init[i]) {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          numsmok[i] = 0;
          cigdyal[i] = 99; //to be calculated later
        } else {
          cigst1[i] = 1;
          endsmoke[i] = 0;
          smokyrs[i] = 0;
          numsmok[i] = 0;
          cigdyal[i] = 0; 
        }
      }
      
      if (cigst1[i-1] == 4) { //current smoker
        if (dice[i] < pr_cess[i]) {
          cigst1[i] = 3;
          endsmoke[i] = 1;
          smokyrs[i] = smokyrs[i-1];
          numsmok[i] = cigdyal[i-1];
          cigdyal[i] = 0;
        } else {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          cigdyal[i] = cigdyal[i-1];
          numsmok[i] = 0;
        }
      }
      
      if ((cigst1[i-1] == 3) && (endsmoke[i-1] < relap_cutoff)){
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (dice[i] < (pr_relap(nrow, endsmoke[i-1] + 1)))  {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          cigdyal[i] = numsmok[i-1];
          numsmok[i] = 0;
        } else {
          cigst1[i] = 3;
          endsmoke[i] = endsmoke[i-1] + 1;
          smokyrs[i] = smokyrs[i-1];
          cigdyal[i] = 0;
          numsmok[i] = numsmok[i-1];
        }
      }
      if ((cigst1[i-1] == 3) && (endsmoke[i-1] >= relap_cutoff)){
        cigst1[i] = 3; 
        endsmoke[i] = endsmoke[i-1] + 1;
        smokyrs[i] = smokyrs[i-1];
        cigdyal[i] = 0;
        numsmok[i] = numsmok[i-1];
      }
    } else { // if different person, start again from year 0
      if ((cigst1[i] == 1) && (dice[i] < pr_init[i])) { //non smoker
        cigst1[i] = 4;
        endsmoke[i] = 0;
        smokyrs[i] = smokyrs[i] + 1;
        numsmok[i] = 0;
        cigdyal[i] = 99; //to be calculated later
      }
      
      if (cigst1[i] == 4) { //current smoker
        if (dice[i] < pr_cess[i]) {
          cigst1[i] = 3;
          endsmoke[i] = 1;
          numsmok[i] = cigdyal[i];
          cigdyal[i] = 0;
        } else {
          smokyrs[i] = smokyrs[i] + 1;
          endsmoke[i] = 0;
          numsmok[i] = 0;
        }
      }
      
      if ((cigst1[i] == 3) && (endsmoke[i] < relap_cutoff)){
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (dice[i] < (pr_relap(nrow, endsmoke[i] + 1)))  {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          cigdyal[i] = numsmok[i];
          smokyrs[i] = smokyrs[i-1] + 1;
          numsmok[i] = 0;
        } 
      }
      if (cigst1[i] == 3) endsmoke[i] = endsmoke[i] + 1;
    }
  }
  return DataFrame::create(_["cigst1"]= cigst1, _["id"]= id, _["year"]= year,
                           _["endsmoke"]= endsmoke, _["smokyrs"]= smokyrs, 
                           _["cigdyal"]= cigdyal, _["numsmok"]= numsmok, 
                           _["cigdyal.rank"]= cigdyal_rank, _["age"]= age,
                           _["sex"]= sex, _["qimd"]= qimd);
}

// [[Rcpp::export]]
DataFrame simsmok(IntegerVector cigst1, NumericVector pr_init, NumericVector pr_cess, NumericVector dice, 
                  IntegerVector id, IntegerVector age, IntegerVector sex, IntegerVector qimd, IntegerVector endsmoke, 
                  IntegerVector smokyrs, NumericVector cigdyal, NumericVector numsmok, IntegerVector year, NumericVector cigdyal_rank,
                  NumericMatrix pr_relap, int relap_cutoff) {
  
  
  // id should be sorted by year
  const int n = cigst1.size();
  int nrow = 0;
  
  if ((cigst1[0] == 1) && (dice[0] < pr_init[0])) { //non smoker
    cigst1[0] = 4;
    endsmoke[0] = 0;
    smokyrs[0] = smokyrs[0] + 1;
    numsmok[0] = 0;
    cigdyal[0] = 99; //to be calculated later
  }
  
  if (cigst1[0] == 4) { //current smoker
    if (dice[0] < pr_cess[0]) {
      cigst1[0] = 3;
      endsmoke[0] = 1;
      numsmok[0] = cigdyal[0];
      cigdyal[0] = 0;
    } else smokyrs[0] = smokyrs[0] + 1; 
  }
  
  if ((cigst1[0] == 3) && (endsmoke[0] < relap_cutoff)){
    switch (sex[0]) {
    case 1: nrow = qimd[0] - 1; break;
    case 2: nrow = qimd[0] + 4; break;
    }
    if (dice[0] < (pr_relap(nrow, endsmoke[0] + 1)))  {
      cigst1[0] = 4;
      endsmoke[0] = 0;
      smokyrs[0] = smokyrs[0] + 1;
      cigdyal[0] = numsmok[0];
      numsmok[0] = 0;
    } 
  }
  if (cigst1[0] == 3) endsmoke[0] = endsmoke[0] + 1;
  
  for (int i = 1; i < n; i++) // start from second element
  {
    if (id[i] == id[i-1]) { // if the same person, take account previous year
      if (cigst1[i-1] == 1) { // never smoker
        if (dice[i] < pr_init[i]) {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          numsmok[i] = 0;
          cigdyal[i] = 99; //to be calculated later
        } else {
          cigst1[i] = 1;
          endsmoke[i] = 0;
          smokyrs[i] = 0;
          numsmok[i] = 0;
          cigdyal[i] = 0; 
        }
      }
      
      if (cigst1[i-1] == 4) { //current smoker
        if (dice[i] < pr_cess[i]) {
          cigst1[i] = 3;
          endsmoke[i] = 1;
          smokyrs[i] = smokyrs[i-1];
          numsmok[i] = cigdyal[i-1];
          cigdyal[i] = 0;
        } else {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          cigdyal[i] = cigdyal[i-1];
          numsmok[i] = 0;
        }
      }
      
      if ((cigst1[i-1] == 3) && (endsmoke[i-1] < relap_cutoff)){
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (dice[i] < (pr_relap(nrow, endsmoke[i-1] + 1)))  {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          smokyrs[i] = smokyrs[i-1] + 1;
          cigdyal[i] = numsmok[i-1];
          numsmok[i] = 0;
        } else {
          cigst1[i] = 3;
          endsmoke[i] = endsmoke[i-1] + 1;
          smokyrs[i] = smokyrs[i-1];
          cigdyal[i] = 0;
          numsmok[i] = numsmok[i-1];
        }
      }
      if ((cigst1[i-1] == 3) && (endsmoke[i-1] >= relap_cutoff)){
        cigst1[i] = 3; 
        endsmoke[i] = endsmoke[i-1] + 1;
        smokyrs[i] = smokyrs[i-1];
        cigdyal[i] = 0;
        numsmok[i] = numsmok[i-1];
      }
    } else { // if different person, start again from year 0
      if ((cigst1[i] == 1) && (dice[i] < pr_init[i])) { //non smoker
        cigst1[i] = 4;
        endsmoke[i] = 0;
        smokyrs[i] = smokyrs[i] + 1;
        numsmok[i] = 0;
        cigdyal[i] = 99; //to be calculated later
      }
      
      if (cigst1[i] == 4) { //current smoker
        if (dice[i] < pr_cess[i]) {
          cigst1[i] = 3;
          endsmoke[i] = 1;
          numsmok[i] = cigdyal[i];
          cigdyal[i] = 0;
        } else {
          smokyrs[i] = smokyrs[i] + 1;
          endsmoke[i] = 0;
          numsmok[i] = 0;
        }
      }
      
      if ((cigst1[i] == 3) && (endsmoke[i] < relap_cutoff)){
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (dice[i] < (pr_relap(nrow, endsmoke[i] + 1)))  {
          cigst1[i] = 4;
          endsmoke[i] = 0;
          cigdyal[i] = numsmok[i];
          smokyrs[i] = smokyrs[i-1] + 1;
          numsmok[i] = 0;
        } 
      }
      if (cigst1[i] == 3) endsmoke[i] = endsmoke[i] + 1;
    }
  }
  return DataFrame::create(_["cigst1"]= cigst1, _["id"]= id, _["year"]= year,
                           _["endsmoke"]= endsmoke, _["smokyrs"]= smokyrs, 
                           _["cigdyal"]= cigdyal, _["numsmok"]= numsmok, 
                           _["cigdyal.rank"]= cigdyal_rank, _["age"]= age,
                           _["sex"]= sex, _["qimd"]= qimd);
}

// [[Rcpp::export]]
NumericVector  numsmok_fix(IntegerVector  cigst1, IntegerVector id, 
                           NumericVector cigdyal, NumericVector numsmok) {
  // id should be sorted and same length as x
  int n = cigst1.size();
  for (int i = 1; i < n; i++) //start loop from 2nd element 
  {
    if (id[i] == id[i-1] && cigst1[i] == 3 && cigst1[i-1] == 4)
    {
      numsmok[i] = cigdyal[i-1];
    } 
    if (id[i] == id[i-1] && cigst1[i] == 3 && cigst1[i-1] == 3)
    {
      numsmok[i] = numsmok[i-1];
    } 
  }
  return numsmok; 
}

// [[Rcpp::export]]
NumericVector shift_byidNum(NumericVector  x, int lag, 
                            float replace, IntegerVector id) {
  // id should be sorted and same length as x
  int n = x.size();
  NumericVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)  
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out; 
  } else {
    stop("This function does not work because number of ids too small!");
  }
}

// [[Rcpp::export]]
IntegerVector shift_byidInt(IntegerVector  x, int lag, 
                            int replace, IntegerVector id) {
  // id should be sorted and same length as x
  int n = x.size();
  IntegerVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)  
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out; 
  } else {
    stop("This function does not work because number of ids too small!");
  }
}

// [[Rcpp::export]]
StringVector  shift_byidStr(CharacterVector x, int lag, 
                            std::string replace, IntegerVector id) {
  // id should be sorted and same length as x
  int n = x.size();
  StringVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)  
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out; 
  } else {
    stop("This function does not work because number of ids too small!");
  }
}

// [[Rcpp::export]]
DataFrame  incidence_events(IntegerVector diab, NumericVector pr_diab_cvd, NumericVector pr_diab_nocvd , 
                            NumericVector dice_diab_inc, IntegerVector chd, NumericVector pr_chd_inc, 
                            NumericVector chd_diab_rr, NumericVector dice_chd_inc, IntegerVector stroke, 
                            NumericVector pr_stroke_inc, NumericVector stroke_diab_rr, NumericVector dice_stroke_inc,
                            IntegerVector id, IntegerVector year, int lag) {
  // id should be sorted and same length as other inputs and year >= 0
  
  const int n = id.size();
  IntegerVector out_diab(n);
  IntegerVector out_chd(n);
  IntegerVector out_stroke(n);
  int i = 0;
  int counter = 0;
  //CHD~diabetes (need to be first for the lag)
  if (diab[i]==2 && year[i] == 0)
  {
    pr_chd_inc[i]    *=  chd_diab_rr[i];
    pr_stroke_inc[i] *=  stroke_diab_rr[i];
  }
  
  //diabetes
  if (diab[i] == 1)
  {  
    if (chd[i] == 0 && stroke[i] == 0)
    {
      out_diab[i] = (dice_diab_inc[i] < pr_diab_nocvd[i]) ? 2 : 1;
    }
    else 
    {
      out_diab[i] = (dice_diab_inc[i] < pr_diab_cvd[i]) ? 2 : 1;
    }
  }
  else out_diab[i] = 2;
  
  // CHD
  if (chd[i] == 0)
  {
    out_chd[i] = (dice_chd_inc[i] < pr_chd_inc[i]) ? 1 : 0;
  }
  if (chd[i] > 0) out_chd[i] = chd[i] + 1;
  //stroke
  if (stroke[i] == 0)
  {
    out_stroke[i] = (dice_stroke_inc[i] < pr_stroke_inc[i]) ? 1 : 0;
  }
  if (stroke[i] > 0) out_stroke[i] = stroke[i] + 1;
  
  
  for (int i = 1; i < n; i++) //start loop from 2nd element 
  { 
    if (id[i] == id[i-1]) // if same person
    {
      counter++ ;//to be used for diabetes lag
      // diabetes
      if (out_diab[i-1] == 1)
      {
        if (out_chd[i-1] == 0 && out_stroke[i-1] == 0) 
        {
          out_diab[i] = (dice_diab_inc[i] < pr_diab_nocvd[i]) ? 2 : 1;
        }
        else
        {
          out_diab[i] = (dice_diab_inc[i] < pr_diab_cvd[i]) ? 2 : 1;
        }
      }
      else out_diab[i] = 2; 
      
      //CVD~diabetes
      if (counter >= lag && out_diab[i-lag] == 2)  
      {
        pr_chd_inc[i]    *=  chd_diab_rr[i];
        pr_stroke_inc[i] *=  stroke_diab_rr[i];
      }
      else if (counter < lag && out_diab[i-counter] == 2)// those diab at year 0 assumed diab for at least as long as cvd.lag
      {
        pr_chd_inc[i]    *=  chd_diab_rr[i];
        pr_stroke_inc[i] *=  stroke_diab_rr[i];
      }
      
      // CHD
      if (out_chd[i-1] == 0)
      {
        out_chd[i] = (dice_chd_inc[i] < pr_chd_inc[i]) ? 1 : 0;
      }
      if (out_chd[i-1] > 0) out_chd[i] = out_chd[i-1] + 1;
      //stroke
      if (out_stroke[i-1] == 0)
      {
        out_stroke[i] = (dice_stroke_inc[i] < pr_stroke_inc[i]) ? 1 : 0;
      }
      if (out_stroke[i-1] > 0) out_stroke[i] = out_stroke[i-1] + 1;
    }
    else // if different person
    {
      counter = 0;
      //CHD~diabetes (need to be first for the lag)
      if (diab[i]==2 && year[i] == 0)
      {
        pr_chd_inc[i]    *=  chd_diab_rr[i];
        pr_stroke_inc[i] *=  stroke_diab_rr[i];
      }
      
      //diabetes
      if (diab[i] == 1)
      {
        if (chd[i] == 0 && stroke[i] == 0)
        {
          out_diab[i] = (dice_diab_inc[i] < pr_diab_nocvd[i]) ? 2 : 1;
        }
        else
        {
          out_diab[i] = (dice_diab_inc[i] < pr_diab_cvd[i]) ? 2 : 1;
        }
      }
      else out_diab[i] = 2;
      
      // CHD
      if (chd[i] == 0)
      {
        out_chd[i] = (dice_chd_inc[i] < pr_chd_inc[i]) ? 1 : 0;
      }
      if (chd[i] > 0) out_chd[i] = chd[i] + 1;
      //stroke
      if (stroke[i] == 0)
      {
        out_stroke[i] = (dice_stroke_inc[i] < pr_stroke_inc[i]) ? 1 : 0;
      }
      if (stroke[i] > 0) out_stroke[i] = stroke[i] + 1;
    }
  }
  
  return DataFrame::create( _["diabtotr"]= out_diab,
                            _["chd.incidence"]= out_chd, 
                            _["stroke.incidence"]= out_stroke);
}

// [[Rcpp::export]]
IntegerVector  mortality_events(NumericVector px_disease,
                                NumericVector px_disease_dice, 
                                IntegerVector id, 
                                int           code) {
  // id should be sorted and same length as other inputs
  const int n = id.size();
  IntegerVector out(n);
  if (px_disease_dice[0] < px_disease[0])
  {
    out[0] = code;
  } 
  else
  {
    out[0] = 0;
  }
  for (int i = 1; i < n; i++) //start loop from 2nd element 
  {
    if ((id[i] == id[i-1]) && (out[i-1] < 0)) 
    {
      out[i] = out[i-1];
    } 
    else
    {
      if (px_disease_dice[i] < px_disease[i])
      {
        out[i] = code;
      } 
      else
      {
        out[i] = 0;
      }
    }
  }
  return out; 
}

// [[Rcpp::export]]
List collect_output(
    IntegerMatrix  x1,// other mort
    IntegerMatrix  x2,// chd mort
    IntegerMatrix  x3,// stroke mort
    IntegerMatrix  x4,// chd incid
    IntegerMatrix  x5,// stroke incid
    IntegerMatrix  x6,// htn incid
    IntegerMatrix  x7,// diab incid
    IntegerMatrix  a0,// age
    NumericMatrix  a1,// qaly by age
    NumericMatrix  a2,// qaly by condition
    IntegerMatrix  a3,// qimd
    NumericMatrix  a4 // cost by condition
) {
  int nrow = x1.nrow(), ncol = x1.ncol();
  IntegerMatrix out(nrow,ncol);
  NumericMatrix qaly(nrow,ncol);
  NumericMatrix cost(nrow,ncol);
  
  int aux[] = {0, 0, 0};
  int age = 0;
  int qimd = 0;
  double aux_qaly = 0.0;
  double aux_cost = 0.0;
  IntegerVector a = IntegerVector::create(-3, -20, -100);
  NumericVector d = NumericVector::create(1, 1, 1);
  
  for (int i = 0; i < nrow; i++) { //loop through rows
    out(i,0) = x1(i,0); //create id column 1
    aux[0] = 0; 
    aux[1] = x1(i,1) + x2(i,1) + x3(i,1); // if dead
    aux[2] = x4(i,1) + x5(i,1) + x6(i,1) + x7(i,1); //if not dead estimate incidence
    // resolve double deaths
    if (aux[1] == -123) aux[1] = Rcpp::sample(a, 1, FALSE, d)[0];
    if (aux[1] == -120) aux[1] = (std::rand() <  0.5) ? -100 : -20; 
    if (aux[1] == -103) aux[1] = (std::rand() <  0.5) ? -100 : -3; 
    if (aux[1] ==  -23) aux[1] = (std::rand() <  0.5) ? -20 : -3; 
    
    // year 0
    out(i,1) = (aux[1] < 0) ? aux[1] : aux[2]; 
    
    //utility
    qaly(i,0) = x1(i,0); //create id column 1
    if (aux[1] == 0) //if alive
    {
      age       = (a0(i,1) <= 100) ? a0(i,1) : 100; // holds the age
      aux_qaly  = a1(age, 1); // holds qaly for age
      aux_qaly  *= (x4(i,1) == 0) ? 1.0 : a2(0,1);
      aux_qaly  *= (x5(i,1) == 0) ? 1.0 : a2(1,1);
      aux_qaly  *= (x6(i,1) == 0) ? 1.0 : a2(2,1);
      aux_qaly  *= (x7(i,1) == 0) ? 1.0 : a2(3,1);
      qaly(i,1) = aux_qaly;
    }
    else qaly(i,1) = 0.0;
    
    //cost
    cost(i,0) = x1(i,0); //create id column 1
    if (aux[1] < 0)
    {
      if (aux[1] == -20)  aux_cost = a4(1,6); //same for all years, qimd
      else if (aux[1] == -3)   aux_cost = a4(1,7); //same for all years, qimd
      else if (aux[1] == -100) aux_cost = 0.0;     //no cost from other death
    }
    else if ((aux[1] == 0) && (aux[2] > 0)) 
    {
      qimd = a3(i,1) - 1;
      aux_cost = 0.0;
      if (x4(i,1) == 400001) aux_cost += a4(qimd,2); 
      if (x4(i,1) == 400000) aux_cost += a4(qimd+5,2);//400000 and 400001 are exclusive
      if (x5(i,1) == 50010) aux_cost += a4(qimd,5); 
      if (x5(i,1) == 50000) aux_cost += a4(qimd+5,5);//50000 and 50010 are exclusive
      if (x6(i,1) == 6000) aux_cost += a4(qimd,4); 
      if (x7(i,1) == 700) aux_cost += a4(qimd,3); 
    }
    else aux_cost = 0.0;
    cost(i,1) = aux_cost;
    
    // loop through columns (concequent years after year 0)
    for (int j = 2; j < ncol; j++) 
    { // matrix first column is id and second year==0   
      aux[0] = out(i,j-1); 
      aux[1] = x1(i,j) + x2(i,j) + x3(i,j); // if 0 = alive
      aux[2] = x4(i,j) + x5(i,j) + x6(i,j) + x7(i,j); //if not dead estimate incidence
      
      // resolve double deaths
      if (aux[1] == -123) aux[1] = Rcpp::sample(a, 1, FALSE, d)[0];
      if (aux[1] == -120) aux[1] = (std::rand() <  0.5) ? -100 : -20; 
      if (aux[1] == -103) aux[1] = (std::rand() <  0.5) ? -100 : -3; 
      if (aux[1] ==  -23) aux[1] = (std::rand() <  0.5) ? -20 : -3; 
      
      if (aux[0] < 0)  out(i,j) =  NA_INTEGER; //those dead remain dead NA after death. NA_integer<x always true
      else out(i,j) = (aux[1] < 0) ? aux[1] : aux[2]; 
      
      //utility
      if (out(i,j) >= 0) //if alive
      {
        age       = (a0(i,j) <= 100) ? a0(i,j) : 100; // holds the age
        aux_qaly  = a1(age, 1); // holds qaly for age
        aux_qaly  *= (x4(i,j) == 0) ? 1.0 : a2(0,1);
        aux_qaly  *= (x5(i,j) == 0) ? 1.0 : a2(1,1);
        aux_qaly  *= (x6(i,j) == 0) ? 1.0 : a2(2,1);
        aux_qaly  *= (x7(i,j) == 0) ? 1.0 : a2(3,1);
        qaly(i,j) = aux_qaly;
      }
      else qaly(i,j) = 0.0;
      
      //cost
      if (out(i,j) != NA_INTEGER && out(i,j) < 0)
      {
        if (aux[1] == -20)  aux_cost = a4(1,6); //same for all years, qimd
        else if (aux[1] == -3)   aux_cost = a4(1,7); //same for all years, qimd
        else if (aux[1] == -100) aux_cost = 0.0;     //no cost from other death
      }
      else if (out(i,j) != NA_INTEGER && out(i,j) > 0)
      {
        qimd = a3(i,j) - 1;
        aux_cost = 0.0;
        if (x4(i,j) == 400001) aux_cost += a4(qimd,2); 
        if (x4(i,j) == 400000) aux_cost += a4(qimd+5,2);//400000 and 400001 are exclusive
        if (x5(i,j) == 50010) aux_cost += a4(qimd,5); 
        if (x5(i,j) == 50000) aux_cost += a4(qimd+5,5);//50000 and 50010 are exclusive
        if (x6(i,j) == 6000) aux_cost += a4(qimd,4); 
        if (x7(i,j) == 700) aux_cost += a4(qimd,3); 
      }
      else aux_cost = 0.0;
      cost(i,j) = aux_cost;
    }
  }
  return Rcpp::List::create(Rcpp::Named("lifecourse", out),
                            Rcpp::Named("utility", qaly),
                            Rcpp::Named("cost", cost));
}

// [[Rcpp::export]]
IntegerVector HC_coverage(IntegerVector id, // Health Checks coverage
                          NumericVector coverage,
                          NumericVector dice
) {
  const int n = id.size();
  IntegerVector invited(n);
  
  invited[0] = (dice[0] < coverage[0]) ? 1 : 0;
  int counter = (invited[0] == 0) ? 0 : 1;
  for (int i = 1; i < n; i++)
  {
    if (id[i] == id[i-1])
    {
      if (counter > 0 && counter < 5)
      {
        counter ++;
        invited[i] = 0;
      }
      else 
      {
        invited[i] = (dice[i] < coverage[i]) ? 1 : 0;
        counter = (invited[i] == 0) ? 0 : 1;
      }
    }
    else //if different person
    {
      invited[i] = (dice[i] < coverage[i]) ? 1 : 0;
      counter = (invited[i] == 0) ? 0 : 1;
    }
  }
  
  return invited;
}

// [[Rcpp::export]]
IntegerVector HC_effect(IntegerVector id, // Health Checks coverage
                        IntegerVector participated,
                        IntegerVector invited
                          
) {
  const int n = id.size();
  IntegerVector HC_effect(n);
  
  HC_effect[0] = (participated[0] == 1) ? 1 : 0;
  for (int i = 1; i < n; i++)
  {
    if (id[i] == id[i-1])
    {
      HC_effect[i] = ((HC_effect[i-1] == 1 && invited[i] == 0) || participated[i] == 1) ? 1 : 0;
    }
    else //if different person
    {
      HC_effect[i] = (participated[i] == 1) ? 1 : 0;
    }
  }
  return HC_effect;
}

// [[Rcpp::export]]
IntegerVector med_effect(IntegerVector id, // Health Checks coverage
                         IntegerVector med_taken,
                         IntegerVector participation_effect
                           
) {
  const int n = id.size();
  IntegerVector med_effect(n);
  
  med_effect[0] = (med_taken[0] == 1) ? 1 : 0;
  for (int i = 1; i < n; i++)
  {
    if (id[i] == id[i-1])
    {
      med_effect[i] = ((med_effect[i-1] == 1 && participation_effect[i] == 1) || med_taken[i] == 1) ? 1 : 0;
    }// In the subsequent HC if person was on statins, remain on statins even if not prescribed
    else //if different person
    {
      med_effect[i] = (med_taken[i] == 1) ? 1 : 0;
    }
  }
  return med_effect;
}

// [[Rcpp::export]]
NumericVector bariatricsurg_effect(IntegerVector id, // Health Checks coverage
                                   NumericVector bmi,
                         IntegerVector behav_effect,
                         IntegerVector participated
                           
) {
  const int n = id.size();
  int flag = -1;
  NumericVector out(n);
  for (int i = 0; i < n; i++)
  {
    if ((participated[i] == 1 && behav_effect[i] == 1 && bmi[i] >= 50.0) || id[i] == flag)
    {
      flag = id[i];
      out[i] = 30.0;
    }
    else 
    {
      out[i] = bmi[i];
    }
  }
  return out;
}
 
 // [[Rcpp::export]]
 IntegerVector smokclinic_effect(IntegerVector id, // Health Checks coverage
                                 IntegerVector cigst1,
                                    IntegerVector behav_effect,
                                    IntegerVector participated
                                      
 ) {
   const int n = id.size();
   int flag = -1;
   IntegerVector out(n);
   for (int i = 0; i < n; i++)
   {
     if ((participated[i] == 1 && behav_effect[i] == 1) || (id[i] == flag && behav_effect[i] == 1))
     {
       flag = id[i];
       out[i] = 3;
     }
     else 
     {
       out[i] = cigst1[i];
     }
   }
   return out;
 }

using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::IntegerVector predict_svylr(mat x, uword e, vec y, rowvec z, vec rank) {
  if (e > 0) x.shed_col(e-1);
  //x.shed_row(e-1);
  x *= y;
  mat cc = repmat(z, x.n_rows, 1);
  cc.each_col() -= x;
  cc.for_each( [](mat::elem_type& val) { val = Rcpp::stats::plogis_0(val, 1, 0); } );  // C++11! lamda function
  umat cc2 = (cc < repmat(rank, 1, cc.n_cols));
  ucolvec out = sum(cc2,1);
  //convert matrix to vector
  Rcpp::IntegerVector tmp = Rcpp::wrap(out);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}
