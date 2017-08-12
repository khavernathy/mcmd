import QtQuick 2.9
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import io.qt.examples.backend 1.0
import QtQuick.Controls.Styles 1.4
import QtQuick.Window 2.2
import QtCharts 2.2

import FileIO 1.0




ApplicationWindow {
    id: root
    visible: true
    width: 1200
    height: 900
    title: qsTr("MCMD :: Monte Carlo & Molecular Dynamics")

    // custom c++ classes
    BackEnd {
        id: backend
        outputLineNumber: 0
    }
    FileIO {
        id: runlogFile
        property int colOffset: 2; // to read data excluding the line # output
        //homeDir: "/Users/douglasfranz" // determined dynamically now
        //source: homeDir+"/mcmd/testzone/runlog.log" // determined dynamically now
        onError: console.log(msg);
    }


    // main visual stuff begins
    SwipeView {
        id: swipeView
        anchors.fill: parent
        currentIndex: tabBar.currentIndex
        onCurrentIndexChanged: {
            //myText.text = myFile.read();
        }

        Page { // 1 :: Input stuff
            Rectangle {
                id: leftcol
                color: "red"
                border.color: "black"
                border.width: 3
                height: parent.height
                width: parent.width/2
                anchors.margins: 5

               Rectangle {
                    width: parent.width
                    height: 50
                    color: "#cdcdcd"
                    Text {
                        font.pixelSize: 30
                        text: "MCMD input parameters"
                    }
                }
/*
                    InputLine {
                        id: firstInpLine
                        y: 50
                        name: "Job name"
                        defaultIn: "test job"
                    }
                    InputLine {
                        y: firstInpLine.y + firstInpLine.height
                        name: "Mode"
                        defaultIn: "mc"
                    }
                    InputLine {
                        y: firstInpLine.y + 2*firstInpLine.height
                        name: "Input atoms"
                        defaultIn: "/Users/douglasfranz/mcmd/testzone/input.pdb"
                    }
                    InputLine {
                        y: firstInpLine.y + 3*firstInpLine.height
                        name: "Potential"
                        defaultIn: "ljes"
                    }
                    InputLine {
                        y: firstInpLine.y + 4*firstInpLine.height
                        name: "Sorbates"
                        defaultIn: "h2_bss"
                    }
                    InputLine {
                        y: firstInpLine.y + 5*firstInpLine.height
                        name: "Fugacity"
                        defaultIn: "h2"
                    }
                    InputLine {
                        y: firstInpLine.y + 6*firstInpLine.height
                        name: "Basis a"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 7*firstInpLine.height
                        name: "Basis b"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 8*firstInpLine.height
                        name: "Basis c"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 9*firstInpLine.height
                        name: "Basis α"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 10*firstInpLine.height
                        name: "Basis β"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 11*firstInpLine.height
                        name: "Basis γ"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 12*firstInpLine.height
                        name: "Ensemble"
                        defaultIn: "uvt"
                    }
                    InputLine {
                        y: firstInpLine.y + 13*firstInpLine.height
                        name: "Corrtime"
                        defaultIn: "1000"
                    }
                    InputLine {
                        y: firstInpLine.y + 14*firstInpLine.height
                        name: "Final step"
                        defaultIn: "10000000"
                    }
                    InputLine {
                        y: firstInpLine.y + 15*firstInpLine.height
                        name: "Temperature"
                        defaultIn: "77.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 16*firstInpLine.height
                        name: "Pressure"
                        defaultIn: "1.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 17*firstInpLine.height
                        name: "Insert factor"
                        defaultIn: "0.667"
                    }
                    InputLine {
                        y: firstInpLine.y + 18*firstInpLine.height
                        name: "Displace factor"
                        defaultIn: "2.5"
                    }
                    InputLine {
                        y: firstInpLine.y + 19*firstInpLine.height
                        name: "Angle rotation factor"
                        defaultIn: "360.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 20*firstInpLine.height
                        name: "Basis γ"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 21*firstInpLine.height
                        name: "Feynman-Hibbs corrections"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 22*firstInpLine.height
                        name: "F-H order"
                        defaultIn: "4"
                    }
                    InputLine {
                        y: firstInpLine.y + 23*firstInpLine.height
                        name: "write_lammps"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 24*firstInpLine.height
                        name: "Auto-reject option"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 25*firstInpLine.height
                        name: "Auto-reject r"
                        defaultIn: "1.78"
                    }
*/
            }
            Rectangle {
                id: midcol
                anchors.left: leftcol.right
                border.color: "black"
                border.width: 3
                color: "blue"
                height: parent.height
                width: parent.width/4
            }
            Rectangle {
                id: rightcol
                anchors.left: midcol.right
                color: "orange"
                border.color: "black"
                border.width: 3
                height: parent.height
                width: parent.width/4
            }

        }
        Page { // 2 :: Runlog output
            id: outputPage

            ScrollView {
                id: runlogOut
                height: parent.height
                width: parent.width
                anchors.fill: parent
                    ScrollBar.vertical: ScrollBar { id: scroller; height: 30; width: parent.width; visible: false;}
                    TextArea {
                        id: outputText
                        color: "white"
                        background: Rectangle { color: "black" }
                        width: 1000
                        height: 780
                        anchors.topMargin: 10
                        font.family: "Monospace"
                        //text: "Heloo world."
                        onTextChanged: {
                            //console.log("the text changed yo")
                            //positionAt: bottom
                        }
                    }

            }

            Component.onCompleted: {
                console.log("scrollview complete");
                //console.log( "WRITE"+ myFile.write("TEST"))
                outputText.text += runlogFile.read();
            }

            Button {
                id: reloadButton
                text: "reload"
                onClicked: {
                    outputText.text += runlogFile.read();
                    //outputText.text += "THE LINE COUNT IS: "+runlogFile.linecount;
                }
                anchors.right: parent.right
                anchors.bottom: parent.bottom
            }

            Timer {
                id: timer
                //property BackEnd backend: BackEnd {}
                interval: 300 // ms
                repeat: true
                running: root.visible
                onTriggered: {
                    //console.time("start data collection");
                    var newText = runlogFile.read();
                    outputText.text += newText;
                    // outputText.contentHeight contains the real height, which grows.
                    //outputText.y = outputText.contentHeight;
                    //scroller.setPosition(outputText.contentHeight)
                    //scroller.position = 0.95; // not quite right, when output gets large it doesn't catch the end
                    scroller.position = (outputText.contentHeight - outputPage.height)/outputText.contentHeight;
                    //console.log("moved scroller, making array");
                    var steps = new Array();
                    var KEs = new Array();
                    var PEs = new Array();
                    var TEs = new Array();
                    var Ds = new Array();
                    var ETemps = new Array();
                    var ITemps = new Array();

                    // loop each NEW line...
                    var allLines = newText.split("\n"); //outputText.text.split("\n");
                    var i=0;
                    while (allLines.length > i) {
                        var thisLine = allLines[i].split(/(\s)/).filter( function(e) { return e.trim().length > 0; } );;
                        if (allLines[i].indexOf("Step") !== -1) {
                            //console.log(allLines[i]);
                            var step = thisLine[1+runlogFile.colOffset];
                            //console.log(step);
                            steps.push(step);
                        }
                        else if (allLines[i].indexOf("KE") !== -1) {
                            var KE = thisLine[2+runlogFile.colOffset];
                            KEs.push(KE);
                        }
                        else if (allLines[i].indexOf("PE") !== -1) {
                            var PE = thisLine[2+runlogFile.colOffset];
                            PEs.push(PE);
                        }
                        else if (allLines[i].indexOf("Total E") !== -1) {
                            var TE = thisLine[3+runlogFile.colOffset];
                            TEs.push(TE);
                        } else if (allLines[i].indexOf("Diffusion c") !== -1) {
                            var Diff = thisLine[5+runlogFile.colOffset];
                            Ds.push(Diff);
                        } else if (allLines[i].indexOf("Emergent T") !== -1) {
                            var ETemp = thisLine[3+runlogFile.colOffset];
                            ETemps.push(ETemp);
                        } else if (allLines[i].indexOf("Instantaneous T") !== -1) {
                            var ITemp = thisLine[3+runlogFile.colOffset];
                            ITemps.push(ITemp);
                        }

                        i++;
                    }
                    //console.log(steps);
                    //console.log(KEs);
                    var laststep = steps[steps.length-1];

                    energychart.updateKE(laststep, KEs[KEs.length -1]);
                    energychart.updatePE(laststep, PEs[PEs.length -1]);
                    energychart.updateTE(laststep, TEs[TEs.length -1]);
                    diffusionchart.updateDiff(laststep, Ds[Ds.length -1]);
                    temperaturechart.updateETemp(laststep, ETemps[ETemps.length -1]);
                    temperaturechart.updateITemp(laststep, ITemps[ITemps.length -1]);

                }


            }
        }

        Page { // 3 :: graph stuff...
            id: graphspage

            Text {
                id: toptitle
                text: "Live data"
                height: 25
                width: parent.width
            }

            ChartView {
                id: energychart
                theme: ChartView.ChartThemeDark
                title: "Energies"

                //anchors.fill: parent
                anchors.left: parent.left
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateKE(x, y) {
                    ke_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }
                function updatePE(x, y) {
                    pe_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }
                function updateTE(x, y) {
                    te_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: ke_line
                    axisX: ValueAxis {
                        min: 0
                        max: 10
                        tickCount: 5
                        titleText: "Step"

                    }

                    axisY: ValueAxis {
                        min: -0.5
                        max: 1.5
                        titleText: "Energy (kJ/mol)"
                    }
                    name: "Kinetic Energy"
                }
                LineSeries {
                    id: pe_line
                    name: "Potential Energy"
                }
                LineSeries {
                    id: te_line
                    name: "Total Energy"
                }
            }

            ChartView {
                id: diffusionchart
                theme: ChartView.ChartThemeBlueCerulean
                title: "Diffusion Coefficient"

                //anchors.fill: parent
                anchors.left: energychart.right
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateDiff(x, y) {
                    diff_line.append(x,y);
                    if (x > diff_line.axisX.max) {
                        diff_line.axisX.max = x;
                    }
                    else if (x < diff_line.axisX.min) {
                        diff_line.axisX.min = x;
                    }
                    if (y > diff_line.axisY.max) {
                        diff_line.axisY.max = y;
                    }
                    else if (y < diff_line.axisY.min) {
                        diff_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: diff_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "D (cm^2 / s)"
                    }
                    name: "D"
                }
            }

            ChartView {
                id: temperaturechart
                theme: ChartView.ChartThemeQt
                title: "Diffusion Coefficient"

                //anchors.fill: parent
                anchors.left: diffusionchart.right
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateETemp(x, y) {
                    temp_line.append(x,y);
                    if (x > temp_line.axisX.max) {
                        temp_line.axisX.max = x;
                    }
                    else if (x < temp_line.axisX.min) {
                        temp_line.axisX.min = x;
                    }
                    if (y > temp_line.axisY.max) {
                        temp_line.axisY.max = y;
                    }
                    else if (y < temp_line.axisY.min) {
                        temp_line.axisY.min = y;
                    }
                }
                function updateITemp(x, y) {
                    itemp_line.append(x,y);
                    if (x > itemp_line.axisX.max) {
                        itemp_line.axisX.max = x;
                    }
                    else if (x < itemp_line.axisX.min) {
                        itemp_line.axisX.min = x;
                    }
                    if (y > itemp_line.axisY.max) {
                        itemp_line.axisY.max = y;
                    }
                    else if (y < itemp_line.axisY.min) {
                        itemp_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: temp_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "Temperature (K)"
                    }
                    name: "Emergent T"
                }
                LineSeries {
                    id: itemp_line
                    name: "Instantaneous T"
                }
            }
        }
        Page {
            Label {
                text: qsTr("Stuff on the 4 page")
            }
        }

        Page {
            Label {
                text: qsTr("more stuff, 5th")
            }
        }
    }

    footer: TabBar {
        id: tabBar
        currentIndex: swipeView.currentIndex
        TabButton {
            text: qsTr("Input setup")
        }
        TabButton {
            text: qsTr("Runlog")
        }
        TabButton {
            text: qsTr("Live graphs")
        }
        TabButton {
            text: qsTr("4th")
        }
        TabButton {
            text: qsTr("5th")
        }
    }

}
