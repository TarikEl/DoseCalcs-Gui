#include "gui/httpdownload.h"
#include "ui_httpdownload.h"
//#include<QtGlobal>

HttpDownload::HttpDownload(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::HttpDownload)
{
    ui->setupUi(this);
    ui->urlEdit->setText("http://qt.com");
    ui->statusLabel->setWordWrap(true);
    ui->downloadButton->setDefault(true);
    ui->quitButton->setAutoDefault(false);

    //progressDialog = new QProgressDialog(this);

    connect(ui->urlEdit, SIGNAL(textChanged(QString)), this, SLOT(enableDownloadButton()));
    //connect(progressDialog, SIGNAL(canceled()), this, SLOT(cancelDownload()));
    //connect(ui->progressBar, SIGNAL(canceled()), this, SLOT(cancelDownload())); // we dont need it because now we dont have the dialog

    ui->progressBar->setValue(0);

}

HttpDownload::~HttpDownload()
{
    delete ui;
}

void HttpDownload::on_downloadButton_clicked() {

//#if QT_NETWORK_LIB
#if QT_NETWORK_LIB

    manager = new QNetworkAccessManager(this);

    // get url
    url = (ui->urlEdit->text()); // we send it from installManagerObj but the user can change it by a new url

    QFileInfo fileInfo(url.path());
    QString fileName = fileInfo.fileName();

    QString filePath = DownloadDir +"/"+fileName;

    if (fileName.isEmpty())
        fileName = "index.html";

    if (QFile::exists(filePath)) { // ask to overwrite the existed file or not
        if (QMessageBox::question(this, tr("HTTP"), tr("There already exists a file called %1 in " "the current directory. Overwrite?").arg(filePath), QMessageBox::Yes|QMessageBox::No, QMessageBox::No) == QMessageBox::No)
            return;
        QFile::remove(filePath);
    }

    file = new QFile(filePath);
    if (!file->open(QIODevice::WriteOnly)) { QMessageBox::information(this, tr("HTTP"), tr("Unable to save the file %1: %2.").arg(filePath).arg(file->errorString()));
        delete file;
        file = 0;
        return;
    }

    // used for progressDialog
    // This will be set true when canceled from progress dialog
    httpRequestAborted = false;

    //progressDialog->setWindowTitle(tr("HTTP"));
    //progressDialog->setLabelText(tr("Downloading %1.").arg(fileName));
    //progressDialog->show();

    // download button disabled after requesting download
    ui->downloadButton->setEnabled(false);

    startRequest(url);
#endif

}

#if QT_NETWORK_LIB
// called from on_downloadButton_clicked() when download button is clicked
void HttpDownload::startRequest(QUrl url)
{
    // get() method posts a request to obtain the contents of the target request and returns a new QNetworkReply object opened for
    // reading which emits the readyRead() signal whenever new data arrives.
    reply = manager->get(QNetworkRequest(url));

    // Whenever more data is received from the network, this readyRead() signal is emitted
    connect(reply, SIGNAL(readyRead()), this, SLOT(httpReadyRead()));

    // Also, downloadProgress() signal is emitted when data is received
    connect(reply, SIGNAL(downloadProgress(qint64,qint64)), this, SLOT(updateDownloadProgress(qint64,qint64)));

    // This signal is emitted when the reply has finished processing.  After this signal is emitted,
    // there will be no more updates to the reply's data or metadata.
    connect(reply, SIGNAL(finished()), this, SLOT(httpDownloadFinished()));
}
#endif

// is a SLOT connected to SIGNAL readyRead() in startRequest() methode
void HttpDownload::httpReadyRead()
{
    // this slot gets called every time the QNetworkReply has new data.
    // We read all of its new data and write it into the file.
    // That way we use less RAM than when reading it at the finished()
    // signal of the QNetworkReply
#if QT_NETWORK_LIB
    if (file)
        file->write(reply->readAll());
#endif
}


// When download finished or canceled, this will be called
// is a SLOT connected to SIGNAL readyRead() in finished() methode
void HttpDownload::httpDownloadFinished()
{
#if QT_NETWORK_LIB
    // when canceled
    if (httpRequestAborted) {
        if (file) {
            file->close();
            file->remove();
            delete file;
            file = 0;
        }
        reply->deleteLater();
        //progressDialog->hide();
        ui->progressBar->setValue(0);
        return;
    }

    // download finished normally
    //progressDialog->hide();
    ui->progressBar->setValue(0);
    file->flush();
    file->close();

    // get redirection url
    QVariant redirectionTarget = reply->attribute(QNetworkRequest::RedirectionTargetAttribute);

    if (reply->error()) {

        file->remove();
        QMessageBox::information(this, tr("HTTP"), tr("Download failed: %1.").arg(reply->errorString()));
        ui->downloadButton->setEnabled(true);

    } else if (!redirectionTarget.isNull()) {
        QUrl newUrl = url.resolved(redirectionTarget.toUrl());
        if (QMessageBox::question(this, tr("HTTP"), tr("Redirect to %1 ?").arg(newUrl.toString()),QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {

            url = newUrl;
            reply->deleteLater();
            file->open(QIODevice::WriteOnly);
            file->resize(0);
            startRequest(url);
            return;
        }
    } else {

        QString fileName = DownloadDir+"/"+QFileInfo(QUrl(url).path()).fileName();
        ui->statusLabel->setText(tr("Downloaded %1 to %2.").arg(fileName).arg(DownloadDir)); // in case of QDir.currentDir()
        ui->downloadButton->setEnabled(true);
    }

    reply->deleteLater();
    reply = 0;
    delete file;
    file = 0;
    manager = 0;
#endif
}

// During the download progress, it can be canceled
// is a SLOT connected to SIGNAL canceled() of progressDialog in constructor() methode
void HttpDownload::cancelDownload()
{
    ui->statusLabel->setText(tr("Download canceled."));
    httpRequestAborted = true;
#if QT_NETWORK_LIB
    reply->abort();
    ui->downloadButton->setEnabled(true);
#endif
}


// called from InstallManagerObj to send the Url of file that we want to download
void HttpDownload::setUrl(QString urll){
#if QT_NETWORK_LIB
    url = urll;
    ui->urlEdit->setText(urll);
#endif
}

// this to set button enabling according to the value of textLine
// is a SLOT connected to SIGNAL textChanged() of urlEdit in constructor() methode
void HttpDownload::enableDownloadButton()
{
    ui->downloadButton->setEnabled(!(ui->urlEdit->text()).isEmpty());
}


// called from InstallManagerObj to send the download directory of a file that we want to download
void HttpDownload::setDownloadDir(QString DownDir){
    DownloadDir = DownDir;
}

// called when the quit is cliked close download dialog
void HttpDownload::on_quitButton_clicked()
{
    this->close();
}


void HttpDownload::on_urlEdit_returnPressed()
{
#if QT_NETWORK_LIB
    on_downloadButton_clicked();
#endif
}

// is a SLOT connected to SIGNAL readyRead() in downloadProgress() methode
void HttpDownload::updateDownloadProgress(qint64 bytesRead, qint64 totalBytes)
{
    if (httpRequestAborted)
        return;

    //progressDialog->setMaximum(totalBytes);
    //progressDialog->setValue(bytesRead);

    ui->progressBar->setMaximum(totalBytes);
    ui->progressBar->setValue(bytesRead);

}


