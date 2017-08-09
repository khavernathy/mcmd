import QtQuick 2.7
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

ApplicationWindow {
    visible: true
    width: 1000
    height: 800
    title: qsTr("MCMD :: Monte Carlo & Molecular Dynamics")

    SwipeView {
        id: swipeView
        anchors.fill: parent
        currentIndex: tabBar.currentIndex

        Page1 {
            Label {
                text: qsTr("Page1 text")
            }
        }

        Page {
            Label {
                text: qsTr("Second page")
                anchors.centerIn: parent

            }
        }
        Page {
            Label {
                text: qsTr("Stuff on the 3rd page")
            }
        }
    }

    footer: TabBar {
        id: tabBar
        currentIndex: swipeView.currentIndex
        TabButton {
            text: qsTr("First")
        }
        TabButton {
            text: qsTr("Second")
        }
        TabButton {
            text: qsTr("3rd")
        }
    }
}
