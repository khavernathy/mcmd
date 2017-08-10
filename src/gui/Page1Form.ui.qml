import QtQuick 2.7
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

Item {
    property alias textField1: textField1
    property alias button1: button1

    RowLayout {
        width: 305
        height: 44
        anchors.horizontalCenterOffset: -52
        anchors.horizontalCenter: parent.horizontalCenter
        anchors.topMargin: 15
        anchors.top: parent.top

        TextField {
            id: textField1
            placeholderText: qsTr("Text Field")
        }

        Button {
            id: button1
            text: qsTr("Press Me")
        }
    }

    MouseArea {
        id: mouseArea
        x: 0
        y: 2
        width: 100
        height: 100
    }

    ProgressBar {
        id: progressBar
        x: 220
        y: 91
        value: 0.5
    }

    Dial {
        id: dial
        x: 456
        y: 2
    }
}
