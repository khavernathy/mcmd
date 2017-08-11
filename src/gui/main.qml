import QtQuick 2.7
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import io.qt.examples.backend 1.0
import QtQuick.Controls.Styles 1.4
import QtQuick.Window 2.2
import FileIO 1.0



ApplicationWindow {
    id: root
    visible: true
    width: 1000
    height: 800
    title: qsTr("MCMD :: Monte Carlo & Molecular Dynamics")

    BackEnd {
        id: backend
        outputLineNumber: 0
    }
    FileIO {
        id: runlogFile
        homeDir: "/Users/douglasfranz"
        source: homeDir+"/mcmd/testzone/runlog.log"
        onError: console.log(msg);
    }
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
            ScrollView {
                id: runlogOut
                anchors.fill: parent
                    TextArea {
                        id: outputText
                        color: "white"
                        background: Rectangle { color: "black" }
                        width: 1000
                        height: 800
                        font.family: "Monospace"
                        //text: "Heloo world."
                        onTextChanged: {
                            console.log("the text changed yo")
                            //positionAt: bottom
                        }
                    }
            }

            Component.onCompleted: {
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
        }

        Page { // 3 :: graph stuff...
            Label {
                text: qsTr("3rd page")
                anchors.centerIn: parent
            }
            TextField {
                text: "Username: "+backend.userName+" ... and number: "+backend.outputLineNumber
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
            text: qsTr("Runlog (output)")
        }
        TabButton {
            text: qsTr("3rd")
        }
        TabButton {
            text: qsTr("4th")
        }
        TabButton {
            text: qsTr("5th")
        }
    }

}
