# Color palettes

get_pals <- function(pal){
    if(!(pal %in% 1:4)){
        stop("No palette matching supplied argument.")
    }
    switch(pal,
           c("#44B3C2", "#F1A94E", "#E45641", "#5D4C46", "#7B8D8E", "#F2EDD8"),
           c("#f4f6af", "#e3c2bd","#d28eca", "#9b7cc3","#656abb"),
           c("#011f4b", "#03396", "#005b96", "#6497b1", "#b3cde0"),
           c("#987543", "#E4C89F", "#BB9866", "#7D5A26", "#57390D",
             "#988443", "#E4D49F", "#BBA766", "#7D6926", "#57460D",
             "#3A376A", "#78769F", "#525082", "#262357", "#13103D",
             "#304A63", "#6A8094", "#466079", "#1D3751", "#0C2339"))
}


