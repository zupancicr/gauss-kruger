<?php 


class Gauss {
                                             
	///Konstante za $racunanje
	private $wgs84_a=6378137.0;	//m
	private $wgs84_a2=40680631590769;
	private $wgs84_b=6356752.314;	//m
	private $wgs84_b2=40408299981544.4;	
	private $wgs84_e2=0.00669438006676466; //e^2
	private $wgs84_e2_=0.00673949681993606; //e'^2;
			
	private $bessel_a=6377397.155;	//m
	private $bessel_a2=40671194472602.1; //a^2
	private $bessel_b=6356078.963;	//m
	private $bessel_b2=40399739783891.2;
	private $bessel_e2=0.00667437217497493;  //e^2
	private $bessel_e2_=0.00671921874158131; //e'^2
	private $bessel_e4=4.45472439300796e-05;
	private $bessel_e6=2.97324885358744e-07;
	private $bessel_e8=1.98445694176601e-09;

	private $dX=-409.520465;
	private $dY=-72.191827;
	private $dZ=-486.872387;
	private $Alfa=1.49625622332431e-05;
	private $Beta=2.65141935723559e-05;
	private $Gama=-5.34282614688910e-05;
	private $dm=-17.919456e-6;

	private $M0 = 0;
	private $M1 = 0;
	private $M2 = 0;

	private $E=4.76916455578838e-12;
	private $D=3.43836164444015e-9;
	private $C=2.64094456224583e-6;
	private $B=0.00252392459157570;
	private $A=1.00503730599692;
	
	private $PI = 3.14159265359;


	function gk2GPS($x, $y, $h)
    {
        $this->M0 = [1.0, sin($this->Gama), -1 * sin($this->Beta)];
        $this->M1 = [-1 * sin($this->Gama), 1, sin($this->Alfa)];
        $this->M2 = [sin($this->Beta), -sin($this->Alfa), 1];


        $y = ($y - 500000) / 0.9999;
        $x = (1 * $x + 5000000) / 0.9999;

        $ab = (1 * $this->bessel_a + 1 * $this->bessel_b);
        $fi0 = (2 * $x) / $ab;

        $dif = 1.0;
        $p1 = $this->bessel_a * (1 - $this->bessel_e2);
        $n = 25;
        while (abs($dif) > 0 && $n > 0) {
            $L = $p1 * ($this->A * $fi0 - $this->B * sin(2 * $fi0) + $this->C * sin(4 * $fi0) - $this->D * sin(6 * $fi0) + $this->E * sin(8 * $fi0));
            $dif = (2 * ($x - $L) / $ab);
            $fi0 = $fi0 + $dif;
            $n--;
        }
        $N = $this->bessel_a / (sqrt(1 - $this->bessel_e2 * pow(sin($fi0), 2)));
        $t = tan($fi0);
        $t2 = pow($t, 2);
        $t4 = pow($t2, 2);
        $cosFi = cos($fi0);
        $ni2 = $this->bessel_e2_ * pow($cosFi, 2);
        $lambda = 0.261799387799149 + ($y / ($N * $cosFi)) - (((1 + 2 * $t2 + $ni2) * pow($y, 3)) / (6 * pow($N, 3) * $cosFi)) + (((5 + 28 * $t2 + 24 * $t4) * pow($y, 5)) / (120 * pow($N, 5) * $cosFi));
        $fi = $fi0 - (($t * (1 + $ni2) * pow($y, 2)) / (2 * pow($N, 2))) + ($t * (5 + 3 * $t2 + 6 * $ni2 - 6 * $ni2 * $t2) * pow($y, 4)) / (24 * pow($N, 4)) - ($t * (61 + 90 * $t2 + 45 * $t4) * pow($y, 6)) / (720 * pow($N, 6));

        $N = $this->bessel_a / (sqrt(1 - $this->bessel_e2 * pow(sin($fi), 2)));
        $X = ($N + $h) * cos($fi) * cos($lambda);
        $Y = ($N + $h) * cos($fi) * sin($lambda);
        $Z = (($this->bessel_b2 / $this->bessel_a2) * $N + $h) * sin($fi);

        $X -= $this->dX;
        $Y -= $this->dY;
        $Z -= $this->dZ;
        $X /= (1 + $this->dm);
        $Y /= (1 + $this->dm);
        $Z /= (1 + $this->dm);

        $X1 = $X - $this->M0[1] * $Y - $this->M0[2] * $Z;
        $Y1 = -1 * $this->M1[0] * $X + $Y - $this->M1[2] * $Z;
        $Z1 = -1 * $this->M2[0] * $X - $this->M2[1] * $Y + $Z;

        $p = sqrt(pow($X1, 2) + pow($Y1, 2));
        $O = atan2($Z1 * $this->wgs84_a, $p * $this->wgs84_b);
        $SinO = sin($O);
        $Sin3O = pow($SinO, 3);
        $CosO = cos($O);
        $Cos3O = pow($CosO, 3);

        $fif = atan2($Z1 + $this->wgs84_e2_ * $this->wgs84_b * $Sin3O, $p - $this->wgs84_e2 * $this->wgs84_a * $Cos3O);
        $lambdaf = atan2($Y1, $X1);

        $N = $this->wgs84_a / sqrt(1 - $this->wgs84_e2 * pow(sin($fif), 2));
        $hf = $p / cos($fif) - $N;

        $fif = ($fif * 180) / $this->PI;
        $lambdaf = ($lambdaf * 180) / $this->PI;

        $retVal = [$fif, $lambdaf, $hf];

        return $retVal;
    }
}
$gauss = new Gauss();
print_r($gauss->gk2GPS(101673.065, 431807.124, 0));
?>

