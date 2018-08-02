#'mGSZ2 GUI
#'
#'mGSZ2 package shiny-based graphical user interface (GUI)
#'\code{startMGSZ2} starts the GUI.
#'
#'@docType methods
#'@name mGSZ2-GUI
#'@rdname mGSZ2-GUI
#'
#'@examples
#'## Start the GUI
#'\dontrun{
#'startMGSZ2();
#'}
#'
#'@importFrom shiny runApp
#'@export startMGSZ2
#'
startMGSZ2 <- function() {
    appDir <- system.file('shiny', 'mGSZ2_GUI', package='mGSZ2');
    if (appDir == '') {
        stop('Could not find GUI directory. Try re-installing `mGSZ2`.',
             call.=FALSE);
    }

    shiny::runApp(appDir);
}
