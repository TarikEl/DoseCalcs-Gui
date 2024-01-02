#ifndef HTTPDOWNLOAD_H
#define HTTPDOWNLOAD_H

#include <QDialog>

#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QMessageBox>
#include <QProgressDialog>

#ifdef QT_NETWORK_LIB

#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QUrl>

#endif

namespace Ui {
class HttpDownload;
}

class HttpDownload : public QDialog
{
    Q_OBJECT

public:

    explicit HttpDownload(QWidget *parent = 0);
    ~HttpDownload();

public:

    void startRequest(QUrl url);
    void setUrl(QString urll);
    void setDownloadDir(QString DownDir);

private slots:

    void on_downloadButton_clicked();

    void on_quitButton_clicked();

    void on_urlEdit_returnPressed();

    // slot for readyRead() signal
    void httpReadyRead();

    // slot for finished() signal from reply
    void httpDownloadFinished();

    // slot for downloadProgress()
    void updateDownloadProgress(qint64, qint64);

    void enableDownloadButton();
    void cancelDownload();

private:

    Ui::HttpDownload *ui;

#ifdef QT_NETWORK_LIB
    QUrl url;
    QNetworkAccessManager *manager;
    QNetworkReply *reply;
#endif

    QProgressDialog *progressDialog;
    QFile *file;
    QString DownloadDir;
    bool httpRequestAborted;
    qint64 fileSize;

};

#endif // HTTPDOWNLOAD_H
