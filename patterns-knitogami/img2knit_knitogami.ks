import cast_ons;
import stockinette;
import parse_knitogami;

with Carrier as 1:{
    img = parse_knitogami.parse_knit_file("patterns/knit_3g_weave.txt");
    //print(img);
    cast_ons.alt_tuck_cast_on(len(img[0]));

    for i in range(0, len(img)):{
        // xfer Loops across to Front Bed;
        for j in range(0, len(img[0])):{
            // if loop is white and on the back bed, bring front
            if img[i][j] == 255:{
                if i > 0:{
                    if img[i][j] != img[i-1][j]:{
                        xfer Back_Needles[j] across to Front Bed;
                    }
                }
            }
            // if loop is black and on front bed, move back
            if img[i][j] == 0:{
                if i > 0:{
                    if img[i][j] != img[i-1][j]:{ //arrays loop around don't forget
                        xfer Front_Needles[j] across to Back Bed;
                    }
                }
                else:{ //don't forget to transfer the first row, assuming it starts on the front bed
                    xfer Front_Needles[j] across to Back Bed;
                }
            }
        }
        in reverse direction:{
            knit Loops;
        }
    }

}